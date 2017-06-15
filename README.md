# Raytracer

A C++ stochastic and bi-directional CPU raytracer (with super-sampling, ambient occlusion, and more)
(Comments are in French)

## Quick Install & Setup ##

On Linux (Tested on ubuntu 16.10, 16.04, 10.04) : 
Go into project directory
Open Terminal

```
>	Make
>	./raytrace
```
The output picture is created in the project directory and named « Image.ppm »

## Output example

![IMG](http://i.imgur.com/N3g8Vvg.png)  
Example with parameters (in file raytracer.ccp) :  
Samples = 4  
GISamples = 4  
SNSamples = 4  
AmbientOcclusionIntensity = 0.5f  
AmbientOcclusion = true  
BiDirectionnal = true  
PhotonsPerLight = 20  
MaxBounces = 3  

## Contributor

Louis A. aka "TheGameIsFixed"  
Contact : thegameisfixed@gmail.com

## License

MIT License

Copyright (c) 2017 TheGameIsFixed

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

