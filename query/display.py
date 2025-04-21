import re
from PIL import Image, ImageDraw, ImageFont

def display(project_name, log_fn):
# Load input
    with open(f"{project_name}_effect.log", "r") as f:
        lines = f.readlines()[1:-2]  # skip timestamp

    # Parse lines
    pattern = re.compile(r"treatment_(\w+)\s+\[(\w+)\]: Causal Effect:\s+([-+]?[0-9]*\.?[0-9]+)")
    data = []

    for line in lines:
        match = pattern.search(line)
        if match:
            treatment, label, effect = match.groups()
            data.append((line.strip(), float(effect), label.lower()))

    # Sort by effect
    data.sort(key=lambda x: x[1], reverse=True)

    # Parameters
    font = ImageFont.load_default()  # or use truetype if you want custom font
    line_height = 20
    padding = 10
    img_width = 400
    img_height = padding * 2 + line_height * len(data)

    # Create image
    img = Image.new("RGB", (img_width, img_height), color="white")
    draw = ImageDraw.Draw(img)

    # Draw each line
    for i, (line, _, label) in enumerate(data):
        y = padding + i * line_height
        color = "red" if label == "pathogenic" else "green"
        draw.text((padding, y), line, fill=color, font=font)

    # Save image
    out_name = f"{project_name}_output.png"
    img.save(out_name)
    log_fn(f"Saved image as {out_name}")
