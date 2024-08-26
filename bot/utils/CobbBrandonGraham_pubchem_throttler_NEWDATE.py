import time
import pubchempy as pcp
import requests

class PubChemThrottler:
    def __init__(self):
        self.base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        self.session = requests.Session()

    def get_throttling_status(self, headers):
        throttling_header = headers.get('X-Throttling-Control', '')
        # Parse the X-Throttling-Control header to extract status information
        # Here we will just extract the service status as an example
        try:
            service_status = throttling_header.split(",")[2].split(":")[1].strip()
            return service_status.lower()
        except (IndexError, ValueError):
            return "unknown"

    def get_compound(self, cid):
        try:
            response = self.session.get(f"{self.base_url}/compound/cid/{cid}/JSON")
            response.raise_for_status()  # Raise an HTTPError on bad status

            # Get throttling status
            status = self.get_throttling_status(response.headers)
            print(f"Service status: {status}")

            if status == "green":
                return response.json()
            elif status == "yellow":
                print("Service is moderately busy. Slowing down requests.")
                time.sleep(1)  # Adjust sleep time as necessary
            elif status == "red":
                print("Service is busy. Slowing down requests significantly.")
                time.sleep(5)  # Adjust sleep time as necessary
            elif status == "black":
                print("Service is overloaded. Waiting longer before retrying.")
                time.sleep(10)  # Adjust sleep time as necessary
            else:
                print("Unable to determine service status. Proceeding with caution.")
                time.sleep(2)

            # Recursive call after adjusting wait time based on status
            return self.get_compound(cid)

        except requests.HTTPError as http_err:
            print(f"HTTP error occurred: {http_err}")
        except Exception as err:
            print(f"An error occurred: {err}")

# Example usage
throttler = PubChemThrottler()
compound_info = throttler.get_compound(2244)  # Example CID
print(compound_info)
