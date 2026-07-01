# Figure Export Workflow

This document describes how to export figures from PowerPoint (A4 slide) to PDF, then split each page into individual PNG files using ImageMagick.

---

## Step 1 — Export PowerPoint to PDF

In PowerPoint, go to **File > Export > Create PDF/XPS Document** and save as PDF. Ensure the slide size is set to **A4 (210 × 297 mm)** before exporting.

---

## Step 2 — Split PDF pages into PNG files

The tool used is [ImageMagick](https://imagemagick.org/). The command reads the PDF at high density (1000 DPI), uses lossless LZW compression, and writes each page as a separate numbered PNG.

### macOS

Install via Homebrew if not already installed:

```bash
brew install imagemagick
```

Run the conversion:

```bash
magick -density 1000 figure_library_survey.pdf -compress lzw +adjoin figure_%02d.png
```

### Linux

Install via your package manager:

```bash
# Debian/Ubuntu
sudo apt install imagemagick

# Fedora/RHEL
sudo dnf install imagemagick
```

Run the conversion (same command as macOS):

```bash
magick -density 1000 figure_library_survey.pdf -compress lzw +adjoin figure_%02d.png
```

> **Note:** On older ImageMagick versions (< 7), replace `magick` with `convert`:
> ```bash
> convert -density 1000 figure_library_survey.pdf -compress lzw +adjoin figure_%02d.png
> ```
>
> PDF processing may be blocked by a security policy. If you get a policy error, edit `/etc/ImageMagick-6/policy.xml` (or `/etc/ImageMagick-7/policy.xml`) and change the PDF rights line from `none` to `read|write`:
> ```xml
> <policy domain="coder" rights="read|write" pattern="PDF" />
> ```

### Windows

Install ImageMagick from [imagemagick.org/script/download.php](https://imagemagick.org/script/download.php) (choose the Windows installer). Ensure **Ghostscript** is also installed from [ghostscript.com](https://www.ghostscript.com/releases/gsdnld.html) — it is required for PDF reading on all platforms.

Run the conversion in Command Prompt or PowerShell from the folder containing the PDF:

```bat
magick -density 1000 figure_library_survey.pdf -compress lzw +adjoin figure_%02d.png
```

> **Tip:** In PowerShell, `%02d` may need escaping. If the output filenames come out wrong, run the command in **Command Prompt** (`cmd.exe`) instead of PowerShell.

---

## Output

Each page of the PDF is written as a zero-padded PNG:

```
figure_00.png   ← page 1
figure_01.png   ← page 2
figure_02.png   ← page 3
...
```

Files are placed in the current working directory. Run the command from the `figure/` directory or adjust the output path as needed.

---

## Parameter reference

| Flag | Effect |
|---|---|
| `-density 1000` | Rasterise at 1000 DPI — suitable for print-quality figures |
| `-compress lzw` | Lossless compression in the output PNG |
| `+adjoin` | Write each page as a separate file instead of a multi-page image |
| `figure_%02d.png` | Output filename pattern; `%02d` gives zero-padded page numbers |
