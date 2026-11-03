"""User Guide: content/user_guide.md, with its images and instruction videos.

A line ``![alt](path/name.png)`` shows static/images/name.png; a line
``@video name.mp4 | Caption`` shows static/videos/name.mp4 under a caption.
"""

from pathlib import Path

import streamlit as st

from src.common.common import page_setup

params = page_setup()
st.title("📖 User Guide")

md_file = Path("content", "user_guide.md")
image_folder = Path("static", "images")
video_folder = Path("static", "videos")


def _flush(block):
    if block:
        st.markdown(''.join(block))
    block.clear()


if not md_file.exists():
    st.error(f"🚨 Error: Could not find {md_file}")
else:
    block = []
    for line in md_file.read_text(encoding="utf-8").splitlines(keepends=True):
        stripped = line.strip()
        if stripped.startswith("@video "):
            _flush(block)
            name, _, caption = stripped[len("@video "):].partition("|")
            path = video_folder / name.strip()
            if path.exists():
                with st.container(border=True):
                    if caption.strip():
                        st.markdown(f"🎬 **{caption.strip()}**")
                    st.video(str(path))
            else:
                st.warning(f"⚠️ Missing video: {name.strip()}")
        elif stripped.startswith("!["):
            _flush(block)
            image_name = stripped[stripped.find("(") + 1:stripped.find(")")].split("/")[-1]
            image_path = image_folder / image_name
            if image_path.exists():
                st.image(str(image_path), width=800)
            else:
                st.warning(f"⚠️ Missing image: {image_name}")
        else:
            block.append(line)
    _flush(block)
