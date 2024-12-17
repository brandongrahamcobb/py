from io import BytesIO
from PIL import Image

import add_watermark
import adjust_hue_and_saturation

def combine(images: list, names: list) -> BytesIO:
    combined_images = []
    for index, (bytes_io, name) in enumerate(zip(images, names)):
        img = Image.open(bytes_io)
        inverted_image = Image.eval(img, lambda x: 255 - x)
        img_bytes = BytesIO()
        inverted_image.save(img_bytes, format='PNG')
        img_bytes_final = add_watermark(img_bytes, watermark_text=name)
        img_final = Image.open(img_bytes_final)
        combined_images.append(img_final)
    widths, heights = zip(*(img.size for img in combined_images))
    total_width = sum(widths)
    max_height = max(heights)
    combined_img = Image.new('RGB', (total_width, max_height))
    x_offset = 0
    for img in combined_images:
        combined_img.paste(img, (x_offset, 0))
        x_offset += img.width
    output = adjust_hue_and_saturation(combined_img, hue_shift=-180, saturation_shift=160)
    output.seek(0)
    return output
