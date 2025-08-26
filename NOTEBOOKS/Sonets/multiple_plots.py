from PIL import Image
import os

# Configuration - using your exact settings
WORKING_DIR = '/home/ignatsonets/AB_plots_v7_big_font_newgrid/'
FILE_LIST = [
    "Age_group_v7_MAGs_combined_plot_big_font_newgrid.png",
 "Stage_of_severity_of_periodontitis_v7_MAGs_combined_plot_big_font_newgrid.png",
    "Rheumatoid_arthritis_v7_MAGs_combined_plot_big_font_newgrid.png",
    "Treatment_of_periodontitis_v7_MAGs_combined_plot_big_font_newgrid.png",
    "BMI_group_v7_MAGs_combined_plot_big_font_newgrid.png",
    "Arterial_hypertension_v7_MAGs_combined_plot_big_font_newgrid.png"
]
GRID_COLS = 3
GRID_ROWS = 2
OUTPUT_FILE = "MAGs_v7_suppl_combined_plots_big_font_newgrid_v3.png"

def main():
    # Change to working directory
    os.chdir(WORKING_DIR)
    print(f"Working in: {os.getcwd()}")
    
    # Verify all files exist
    missing = [f for f in FILE_LIST if not os.path.exists(f)]
    if missing:
        print(f"Error: Missing files - {', '.join(missing)}")
        return
    
    # Load images and get dimensions from first file
    images = []
    for f in FILE_LIST:
        img = Image.open(f)
        images.append(img)
    
    # Get dimensions from first image
    width, height = images[0].size
    print(f"Detected dimensions: {width}×{height} pixels")
    
    # Create blank canvas
    combined = Image.new('RGB', (GRID_COLS * width, GRID_ROWS * height))
    
    # Arrange in grid (left to right, top to bottom)
    for i, img in enumerate(images):
        row = i // GRID_COLS  # 0-1
        col = i % GRID_COLS    # 0-2
        x = col * width
        y = row * height
        combined.paste(img, (x, y))
        print(f"Placed {FILE_LIST[i]} at position ({row},{col})")
    
    # Save result
    combined.save(OUTPUT_FILE)
    print(f"\nSuccess! Combined image saved to:\n{os.path.join(WORKING_DIR, OUTPUT_FILE)}")
    print(f"Final dimensions: {combined.width}×{combined.height} pixels")

if __name__ == "__main__":
    main()
