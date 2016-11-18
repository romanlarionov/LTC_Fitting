
Linearly Transformed Cosines Offline Fitting Code
================================================

This is the code base for the BRDF alignment step of the [LTC paper](https://drive.google.com/file/d/0BzvWIdpUpRx_d09ndGVjNVJzZjA/view).

A link to the main research website which contains the original code as well as a javascript demo: https://eheitzresearch.wordpress.com/415-2/

The code in this repository does not drastically change the fitting process. It does include additional code for JSON and PNG formats.

###One small note

If you're having trouble compiling with the given makefile, you might need to download and install [ImageMagick](http://www.imagemagick.org/script/index.php) and [GraphicsMagick](http://www.graphicsmagick.org/). A simple brew install worked on my mac.

Full credit goes to the team at Unity for writing this code.

MERL Database
============

If you want the additional code for fitting an LTC to a [MERL](http://www.merl.com/brdf/) material type, refer to the `merl` branch of this repository.
