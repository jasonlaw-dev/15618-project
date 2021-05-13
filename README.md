# Canny Edge Detection in OpenMPI

## How to run

Note that this is only tested on the GHC machines.

1. Download [FreeImage](http://downloads.sourceforge.net/freeimage/FreeImage3180.zip)

2. Place the zip file **next to** the project directory. E.g.
  ```
  dir/
    15618-project
    FreeImage3180.zip
  ```

3. Build FreeImage
  ```sh
  unzip FreeImage3180.zip
  cd FreeImage
  make
  ```

4. Build the project
  ```sh
  cd ../15618-project
  make
  ```
5. Run the program using mpirun
  ```sh
  mpirun -np 1 ./main -f images/wallpaper.jpg # this runs the bilateral filter pipeline with 1 core on wallpaper.jpg
  ```
  Flags supported:
   - `-f [image]` Specify the image file to be generated
   - `-s [stage]` Specify running the pipeline up to a specific stage, accepts 1-5. This only prints the runtime for that particular stage. When not specified, it will run to the end and print the overall runtime.
   - `-g` Change the blur algorithm from bilateral filter to Gaussian blur
   - `-t` Set the Gaussian blur to two-pass, only works when `-g` is used
   - `-o` Turn on optimization mode

6. This generates an image by appending the stage name to the filename, e.g. `wallpaper-5-hyster.jpg`

## List of images
- images/audrey.png
- images/bigben.png
- images/chess.png - the smallest image
- images/emma.png
- images/panther.jpg - the largest image
- images/valve.png
- images/wallpaper.jpg - our benchmark input