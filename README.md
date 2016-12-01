
Linearly Transformed Cosines Offline Fitting Code
================================================

This is the code base for the BRDF alignment step of the [LTC paper](https://drive.google.com/file/d/0BzvWIdpUpRx_d09ndGVjNVJzZjA/view).

A link to the main research website which contains the original code as well as a javascript demo: https://eheitzresearch.wordpress.com/415-2/

The code in this repository does not drastically change the fitting process. It does include code for additional export procedures and BRDFs to fit. The current implimentation replaces the importance sampling process with a uniform sampling approach. This is primarily due to the fact that constructing a PDF/CDF for the MERL 2d data would require more advanced algorithms than I am currently prepared to deal with. It should also be noted that the type of MERL data should be carefully chosen to be isotropic. Anisotropic BRDFs are not supported by this fitting code. A safe bet would be any metallic material.

###One small note

If you're having trouble compiling with the given makefile, you might need to download and install [ImageMagick](http://www.imagemagick.org/script/index.php) and [GraphicsMagick](http://www.graphicsmagick.org/). A simple brew install worked on my mac.

Full credit goes to the team at Unity for writing this code.
