def gsrs(arg):
    chrome_options = Options()
    chrome_options.add_argument('--headless')  # Run headless Chrome (no UI)
    chrome_options.add_argument('--no-sandbox')
    chrome_options.add_argument('--disable-dev-shm-usage')
    driver = webdriver.Chrome(service=Service(executable_path='chromedriver'), options=chrome_options)
    try:
        search_url = f'https://gsrs.ncats.nih.gov/ginas/app/beta/browse-substance?search={arg}'
        driver.get(search_url)
        driver.implicitly_wait(10)  # Adjust the wait time as needed
        img_element = driver.find_element(By.CSS_SELECTOR, 'body > app-root > app-base > app-substances-browse > div > div.substance-cards > app-substance-summary-card > mat-card > mat-card-title > a')
        if img_element:
            img_src = img_element.get_attribute('href')
            if img_src:
                stripped = img_src.split('/', -1)[-1:]
                link = f'https://gsrs.ncats.nih.gov/api/v1/substances/render({stripped[0]})?format=png&size=512&stereo=true'
                response = requests.get(link)
                image_bytes = response.content
                image = Image.open(BytesIO(image_bytes))
                image = image.convert('RGBA')
                draw = ImageDraw.Draw(image)
                width, height = image.size
                diagonal = math.sqrt(width**2 + height**2)
                font_size = int(diagonal / 15)
                try:
                    font = ImageFont.truetype('Roboto-Regular.ttf', font_size)
                except IOError:
                    font = ImageFont.load_default()
                bbox = draw.textbbox((0, 0), arg, font=font)
                text_width = bbox[2] - bbox[0]
                text_height = bbox[3] - bbox[1]
                text_x = (width - text_width) / 2
                text_y = (height - text_height) / 2
                watermark_image = Image.new('RGBA', image.size, (0, 0, 0, 0))
                watermark_draw = ImageDraw.Draw(watermark_image)
                watermark_draw.text((text_x, text_y), arg, font=font, fill=(255, 255, 255, 64))
                mask = watermark_image.split()[3]
                image.paste(watermark_image, (0, 0), mask)
                return image
            else:
                return 'No src attribute found in the <img> element'
        else:
            return 'No <img> element found with the specified CSS path'
    finally:
        driver.quit()
