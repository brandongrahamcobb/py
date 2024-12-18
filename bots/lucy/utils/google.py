from bs4 import BeautifulSoup

import requests

def google(query: str, num_results: int = 5):
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
    }
    search_url = "https://www.google.com/search"
    params = {"q": query, "num": num_results}
    try:
        response = requests.get(search_url, headers=headers, params=params)
        response.raise_for_status()  # Raise an error for bad HTTP responses
        soup = BeautifulSoup(response.text, "html.parser")
        results = []
        for g in soup.find_all("div", class_="g"):  # Updated class name
            title = g.find("h3").text if g.find("h3") else "No title"
            link = g.find("a")["href"] if g.find("a") else "No link"
            results.append({"title": title, "link": link})
            if len(results) >= num_results:
                break
        return results
    except requests.exceptions.RequestException as e:
        print(f"Error during the web request: {e}")
        return []

