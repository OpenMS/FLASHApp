import streamlit as st
from pathlib import Path
from src.common.common import page_setup

# Setup page
params = page_setup()
st.title("📖 User Guide")

# Define paths
md_file = Path("content", "user_guide.md")
image_folder = Path("static", "Images")

# Read the User Guide Markdown file
if md_file.exists():
    with open(md_file, "r", encoding="utf-8") as f:
        content = f.readlines()
else:
    st.error(f"🚨 Error: Could not find {md_file}")
    content = []

# Process Markdown content and replace image references
for line in content:
    # Custom image placeholder for flashdeconv_configure.png
    if "[FLASHDECONV_CONFIGURE_IMAGE]" in line:
        image_path = image_folder / "flashdeconv_configure.png"
        if image_path.exists():
            st.image(str(image_path), caption="FLASHDeconv Configuration", width=300)
        else:
            st.warning("⚠️ Image 'flashdeconv_configure.png' not found.")

    # Markdown-style image (e.g., ![alt text](path/to/image.png))
    elif line.strip().startswith("!["):
        start = line.find("(") + 1
        end = line.find(")")
        image_name = line[start:end].split("/")[-1]
        image_path = image_folder / image_name

        if image_path.exists():
            st.image(str(image_path), caption=image_name, width=800)
        else:
            st.warning(f"⚠️ Missing image: {image_name}")

    # Normal markdown text
    else:
        st.markdown(line)
