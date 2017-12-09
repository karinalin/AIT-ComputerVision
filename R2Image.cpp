// Source file for image class


// Include files 

#include "R2/R2.h"
#include "R2Pixel.h"
#include "R2Image.h"
#include "svd.h"



////////////////////////////////////////////////////////////////////////
// Constructors/Destructors
////////////////////////////////////////////////////////////////////////

R2Image::
R2Image(void)
  : pixels(NULL),
    npixels(0),
    width(0), 
    height(0)
{
}



R2Image::
R2Image(const char *filename)
  : pixels(NULL),
    npixels(0),
    width(0), 
    height(0)
{
  // Read image
  Read(filename);
}



R2Image::
R2Image(int width, int height)
  : pixels(NULL),
    npixels(width * height),
    width(width), 
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);
}



R2Image::
R2Image(int width, int height, const R2Pixel *p)
  : pixels(NULL),
    npixels(width * height),
    width(width), 
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = p[i];
}



R2Image::
R2Image(const R2Image& image)
  : pixels(NULL),
    npixels(image.npixels),
    width(image.width), 
    height(image.height)
    
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = image.pixels[i];
}



R2Image::
~R2Image(void)
{
  // Free image pixels
  if (pixels) delete [] pixels;
}



R2Image& R2Image::
operator=(const R2Image& image)
{
  // Delete previous pixels
  if (pixels) { delete [] pixels; pixels = NULL; }

  // Reset width and height
  npixels = image.npixels;
  width = image.width;
  height = image.height;

  // Allocate new pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = image.pixels[i];

  // Return image
  return *this;
}

struct Feature
{
    int centerX;
    int centerY;
    int endX;
    int endY;
    R2Pixel HarrisValue;
    
    Feature(int x, int y, R2Pixel val)
    {
        centerX = x;
        centerY = y;
        endX = 0;
        endY = 0;
        HarrisValue = val;
    }
    
    bool operator<(const Feature& feature) const {
        double valueIntensity = HarrisValue[0] + HarrisValue[1] + HarrisValue[2];
        double featureIntensity = feature.HarrisValue[0] + feature.HarrisValue[1] + feature.HarrisValue[2];
        
        return valueIntensity < featureIntensity;
    }
    
};

double** R2Image::
svdTest(std::vector<double> svdTestInput)
{

  //PROBLEM 1
  //Store pairs of points as feature vectors
  std::vector<Feature> pointList;

  int inputLength = svdTestInput.size()/4;
  //fprintf(stderr, "inputLength: %d\n", inputLength);

  for (int i=0; i<inputLength; i++) {

    double x1 = svdTestInput[(i*4)];
    double y1 = svdTestInput[(i*4) + 1];
    double x2 = svdTestInput[(i*4) + 2];
    double y2 = svdTestInput[(i*4) + 3];

    //fprintf(stderr, "x1: %f, y1: %f, x2: %f, y2: %f\n", x1, y1, x2, y2);

    Feature tempPoint = Feature(x1, y1, R2black_pixel);
    tempPoint.endX = x2;
    tempPoint.endY = y2;
    pointList.push_back(tempPoint);
  }

  //Build the 9x8 matrix of equations
  double** decomposeMatrix = dmatrix(1, inputLength*2, 1, 9);
  double** submatrix = dmatrix(1,2,1,9);

  submatrix[1][1] = 0;
  submatrix[1][2] = 0;
  submatrix[1][3] = 0;
  submatrix[1][6] = -1;
  submatrix[2][3] = 1;
  submatrix[2][4] = 0;
  submatrix[2][5] = 0;
  submatrix[2][6] = 0;

  //fprintf(stderr, "%s\n", "static values added to submatrix");

  for (int i = 0; i < inputLength; i++) {
    //Build 9x2 submatrix for each point pair
    double x1 = pointList[i].centerX;
    double y1 = pointList[i].centerY;
    double x2 = pointList[i].endX;
    double y2 = pointList[i].endY;

    //fprintf(stderr, "x1: %f, y1: %f, x2: %f, y2: %f\n", x1, y1, x2, y2);

    submatrix[1][4] = -1 * x1;
    submatrix[1][5] = -1 * y1;
    submatrix[1][7] = x1*y2;
    submatrix[1][8] = y1*y2;
    submatrix[1][9] = y2;
    submatrix[2][1] = x1;
    submatrix[2][2] = y1;
    submatrix[2][7] = -1 * x1*y2;
    submatrix[2][8] = -1 * y1*x2;
    submatrix[2][9] = -1 * x2;

    //fprintf(stderr, "%s\n", "dynamic values added to submatrix");


    //Add submatrices to main matrix to decompose
    for (int b = 1; b < 3; b++) {
      for (int a = 1; a < 10; a++) {
        //fprintf(stderr, "%s\n", "transfering to decomposeMatrix");

        decomposeMatrix[i + i + b][a] = submatrix[b][a];
      }
    }
    printf("sub matrix: \n[%f %f %f %f %f %f %f %f %f \n %f %f %f %f %f %f %f %f %f]\n", submatrix[1][1], submatrix[1][2], submatrix[1][3], submatrix[1][4], submatrix[1][5], submatrix[1][6], submatrix[1][7], submatrix[1][8], submatrix[1][9],  submatrix[2][1], submatrix[2][2], submatrix[2][3], submatrix[2][4], submatrix[2][5], submatrix[2][6], submatrix[2][7], submatrix[2][8], submatrix[2][9]);

  }

    printf("decomposed matrix: \n[%f %f %f %f %f %f %f %f %f \n %f %f %f %f %f %f %f %f %f \n %f %f %f %f %f %f %f %f %f \n %f %f %f %f %f %f %f %f %f \n %f %f %f %f %f %f %f %f %f \n %f %f %f %f %f %f %f %f %f \n %f %f %f %f %f %f %f %f %f \n %f %f %f %f %f %f %f %f %f]\n", 
        decomposeMatrix[1][1], decomposeMatrix[1][2], decomposeMatrix[1][3], decomposeMatrix[1][4], decomposeMatrix[1][5], decomposeMatrix[1][6], decomposeMatrix[1][7], decomposeMatrix[1][8], decomposeMatrix[1][9],  
        decomposeMatrix[2][1], decomposeMatrix[2][2], decomposeMatrix[2][3], decomposeMatrix[2][4], decomposeMatrix[2][5], decomposeMatrix[2][6], decomposeMatrix[2][7], decomposeMatrix[2][8], decomposeMatrix[2][9],

        decomposeMatrix[3][1], decomposeMatrix[3][2], decomposeMatrix[3][3], decomposeMatrix[3][4], decomposeMatrix[3][5], decomposeMatrix[3][6], decomposeMatrix[3][7], decomposeMatrix[3][8], decomposeMatrix[3][9],  
        decomposeMatrix[4][1], decomposeMatrix[4][2], decomposeMatrix[4][3], decomposeMatrix[4][4], decomposeMatrix[4][5], decomposeMatrix[4][6], decomposeMatrix[4][7], decomposeMatrix[4][8], decomposeMatrix[4][9],
        
        decomposeMatrix[5][1], decomposeMatrix[5][2], decomposeMatrix[5][3], decomposeMatrix[5][4], decomposeMatrix[5][5], decomposeMatrix[5][6], decomposeMatrix[5][7], decomposeMatrix[5][8], decomposeMatrix[5][9],  
        decomposeMatrix[6][1], decomposeMatrix[6][2], decomposeMatrix[6][3], decomposeMatrix[6][4], decomposeMatrix[6][5], decomposeMatrix[6][6], decomposeMatrix[6][7], decomposeMatrix[6][8], decomposeMatrix[6][9],
        
        decomposeMatrix[7][1], decomposeMatrix[7][2], decomposeMatrix[7][3], decomposeMatrix[7][4], decomposeMatrix[7][5], decomposeMatrix[7][6], decomposeMatrix[7][7], decomposeMatrix[7][8], decomposeMatrix[7][9],  
        decomposeMatrix[8][1], decomposeMatrix[8][2], decomposeMatrix[8][3], decomposeMatrix[8][4], decomposeMatrix[8][5], decomposeMatrix[8][6], decomposeMatrix[8][7], decomposeMatrix[8][8], decomposeMatrix[8][9]);


  printf("\n PROBLEM 1:\n");


  // Compute the SVD
  double singularValues[10];
  double** nullspaceMatrix = dmatrix(1,9,1,9);
  svdcmp(decomposeMatrix, inputLength * 2, 9, singularValues, nullspaceMatrix);

  // Get the result
  printf("\n Singular values: %f, %f, %f, %f, %f, %f, %f, %f, %f\n", singularValues[1], singularValues[2], singularValues[3], singularValues[4], singularValues[5], singularValues[6], singularValues[7], singularValues[8], singularValues[9]);

  // Find the smallest singular value:
  int smallestIndex = 1;
  for(int i=2;i<10;i++) if(singularValues[i]<singularValues[smallestIndex]) smallestIndex=i;

  // Solution for H matrix is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
  printf("Conic coefficients: %f, %f, %f, %f, %f, %f, %f, %f, %f\n", nullspaceMatrix[1][smallestIndex], nullspaceMatrix[2][smallestIndex], nullspaceMatrix[3][smallestIndex], nullspaceMatrix[4][smallestIndex], nullspaceMatrix[5][smallestIndex], nullspaceMatrix[6][smallestIndex],
    nullspaceMatrix[7][smallestIndex], nullspaceMatrix[8][smallestIndex], nullspaceMatrix[9][smallestIndex]);

  //Compute scalar for H matrix
  //double scalar = 1 / nullspaceMatrix[9][smallestIndex];
  //Scale H matrix by multiplying by above scalar
  //for (int i = 1; i < 10; i++) {
    //nullspaceMatrix[i][smallestIndex] = nullspaceMatrix[i][smallestIndex] * scalar;
  //}

  // Turn solution into 3x3 matrix as final solution
  double** hMatrix = dmatrix(1, 3, 1, 3);
  hMatrix[1][1] = (-1)* nullspaceMatrix[1][smallestIndex];
  hMatrix[1][2] = (-1)* nullspaceMatrix[2][smallestIndex];
  hMatrix[1][3] = (-1)* nullspaceMatrix[3][smallestIndex];
  hMatrix[2][1] = (-1)* nullspaceMatrix[4][smallestIndex];
  hMatrix[2][2] = (-1)* nullspaceMatrix[5][smallestIndex];
  hMatrix[2][3] = (-1)* nullspaceMatrix[6][smallestIndex];
  hMatrix[3][1] = (-1)* nullspaceMatrix[7][smallestIndex];
  hMatrix[3][2] = (-1)* nullspaceMatrix[8][smallestIndex];
  hMatrix[3][3] = (-1)* nullspaceMatrix[9][smallestIndex];

  //Print final solution for H matrix
  printf("H matrix: \n[%f %f %f\n %f %f %f\n %f %f %f]\n", hMatrix[1][1], hMatrix[1][2], hMatrix[1][3], hMatrix[2][1], hMatrix[2][2], hMatrix[2][3], hMatrix[3][1], hMatrix[3][2], hMatrix[3][3]);

  return hMatrix;
}



////////////////////////////////////////////////////////////////////////
// Image processing functions
// YOU IMPLEMENT THE FUNCTIONS IN THIS SECTION
////////////////////////////////////////////////////////////////////////

// Per-pixel Operations ////////////////////////////////////////////////

void R2Image::
Brighten(double factor)
{
  // Brighten the image by multiplying each pixel component by the factor.
  // This is implemented for you as an example of how to access and set pixels
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      Pixel(i,j) *= factor;
      Pixel(i,j).Clamp();
    }
  }
}

void R2Image::
SobelX(void)
{
  // Apply the Sobel oprator to the image in X direction
  
  R2Image* tempImage = new R2Image(width, height);

  for (int i = 1; i < width-1; i++) {
    for (int j = 1; j < height-1; j++) {

      tempImage->Pixel(i,j) += Pixel(i-1, j+1) * 1.0;
      tempImage->Pixel(i,j) += Pixel(i-1,j) * 2.0;
      tempImage->Pixel(i,j) += Pixel(i-1,j-1) * 1.0;
      tempImage->Pixel(i,j) += Pixel(i+1,j+1) * -1.0;
      tempImage->Pixel(i,j) += Pixel(i+1,j) * -2.0;
      tempImage->Pixel(i,j) += Pixel(i+1,j-1) * -1.0;

      tempImage->Pixel(i,j);//.Clamp();
    }
  }

  for (int k = 1; k < width-1; k++) {
    for (int p = 1; p < height-1; p++) {

      Pixel(k,p) = tempImage ->Pixel(k,p);
    }
  }
}

void R2Image::
SobelY(void)
{
  // Apply the Sobel oprator to the image in Y direction
  
  R2Image* tempImage = new R2Image(width, height);

  for (int i = 1; i < width-1; i++) {
    for (int j = 1; j < height-1; j++) {

      tempImage->Pixel(i,j) += Pixel(i-1, j+1) * -1.0;
      tempImage->Pixel(i,j) += Pixel(i,j+1) * -2.0;
      tempImage->Pixel(i,j) += Pixel(i+1,j+1) * -1.0;
      tempImage->Pixel(i,j) += Pixel(i-1,j-1) * 1.0;
      tempImage->Pixel(i,j) += Pixel(i,j-1) * 2.0;
      tempImage->Pixel(i,j) += Pixel(i+1,j-1) * 1.0;

      tempImage->Pixel(i,j);//.Clamp();

    }
  }

  for (int k = 1; k < width-1; k++) {
    for (int p = 1; p < height-1; p++) {

      Pixel(k,p) = tempImage ->Pixel(k,p);
    }
  }
}

void R2Image::
LoG(void)
{
  // Apply the LoG oprator to the image
  
  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  fprintf(stderr, "LoG() not implemented\n");
}



void R2Image::
ChangeSaturation(double factor)
{
  // Changes the saturation of an image
  // Find a formula that changes the saturation without affecting the image brightness

  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  fprintf(stderr, "ChangeSaturation(%g) not implemented\n", factor);
}


// Linear filtering ////////////////////////////////////////////////
void R2Image::
Blur(double sigma)
{
  // Gaussian blur of the image. Separable solution is preferred
  
    //fprintf(stderr, "STARTING BLUR\n");

R2Image* tempImage = new R2Image(width, height);

  R2Pixel outputPixel;

  int kernelLength = 6 * (int)sigma + 1;
  double * kernel = (double*)malloc(kernelLength*sizeof(float));  
  double kernelSum = 0;

  //fprintf(stderr, "looping to do gaussian calculation\n");

  for (int k = 0; k <= kernelLength/2; k++) {
    double gaussianCalculation = ((1)/(sigma * sqrt(2*3.1415926535)))*exp((-pow((-kernelLength/2+k),2))/(2*(pow(sigma, 2))));
    //printf("%d", -kernelLength/2+k);

    kernel[k] = gaussianCalculation;
    kernel[kernelLength-1-k] = gaussianCalculation;
    
  }

  if (kernelLength%2 == 1) {

    double centerCalculation = ((1)/(sigma * sqrt(2*3.1415926535)))*exp((-pow((0),2))/(2*(pow(sigma, 2))));
    kernel[kernelLength/2] = centerCalculation;
  }

  
  //fprintf(stderr, "adding to KernelSum\n");

  for (int y=0; y<kernelLength; y++) {
    kernelSum += kernel[y];

  }

  /*for (int x=0; x<kernelLength; x++) {
    kernel[x] = kernel[x]/kernelSum;
  }*/

  // for (int m=0; m<kernelLength;m++) {

  //   printf("kernel[%d]: %f \n", m, kernel[m]);   
  // }

  //printf("Kernel sum = %f", kernelSum);

  //fprintf(stderr, "adding to output pixel then assigning\n");


  for (int i = kernelLength/2; i < width-(kernelLength/2); i++) {
    for (int j = kernelLength/2; j < height-kernelLength/2; j++) {

      outputPixel.Reset(0.0,0.0,0.0,1.0);

      for (int p = 0; p < kernelLength; p++) {
        outputPixel += Pixel(i-(kernelLength/2)+p, j) * kernel[p];

      }

        tempImage->Pixel(i,j) = outputPixel;
        //tempImage->Pixel(i,j).Clamp();
    }
  }

  //fprintf(stderr, "temp image to original image\n");

  for (int k = kernelLength/2; k < width-(kernelLength/2); k++) {
    for (int p = kernelLength/2; p < height-kernelLength/2; p++) {

      Pixel(k,p) = tempImage ->Pixel(k,p);
    }
  }

 // fprintf(stderr, "...\n");

  for (int q = kernelLength/2; q < width-(kernelLength/2); q++) {
    for (int h = kernelLength/2; h < height-kernelLength/2; h++) {

      outputPixel.Reset(0.0,0.0,0.0,1.0);

      for (int t = 0; t < kernelLength; t++) {
        outputPixel += Pixel(q,h-(kernelLength/2)+t) * kernel[t];
      }
        tempImage->Pixel(q,h) = outputPixel;
        //tempImage->Pixel(q,h).Clamp();
    }
  }

  for (int k = kernelLength/2; k < width-(kernelLength/2); k++) {
    for (int p = kernelLength/2; p < height-kernelLength/2; p++) {

      Pixel(k,p) = tempImage ->Pixel(k,p);
    }
  }
  //free(kernel);
}


void R2Image::line(int x0, int x1, int y0, int y1, float r, float g, float b)
{
  if(x0>x1)
  {
    int x=y1;
    y1=y0;
    y0=x;

    x=x1;
    x1=x0;
    x0=x;
  }
     int deltax = x1 - x0;
     int deltay = y1 - y0;
     float error = 0;
     float deltaerr = 0.0;
   if(deltax!=0) deltaerr =fabs(float(float(deltay) / deltax));    // Assume deltax != 0 (line is not vertical),
           // note that this division needs to be done in a way that preserves the fractional part
     int y = y0;
     for(int x=x0;x<=x1;x++)
   {
     Pixel(x,y).Reset(r,g,b,1.0);
         error = error + deltaerr;
         if(error>=0.5)
     {
       if(deltay>0) y = y + 1;
       else y = y - 1;

             error = error - 1.0;
     }
   }
   if(x0>3 && x0<width-3 && y0>3 && y0<height-3)
   {
     for(int x=x0-3;x<=x0+3;x++)
     {
       for(int y=y0-3;y<=y0+3;y++)
       {
         Pixel(x,y).Reset(r,g,b,1.0);
       }
     }
   }
}



std::vector<Feature> topFeatureVec;
std::vector<Feature> featuresToCheck;
int numFeaturesToCheck;

void R2Image::
Harris(double sigma)
{
    // Harris corner detector. Make use of the previously developed filters, such as the Gaussian blur filter
  // Output should be 50% grey at flat regions, white at corners and black/dark near edges

  //find 150 features with very high corner score
  //make sure there is at least 10 pixel distance between them
  //fprintf(stderr, "height*width: %d\n" , height*width);
  fprintf(stderr, "copying this into 4 images\n");

  R2Image sobelX2(*this); 
  R2Image sobelY2(*this);
  R2Image sobelXY(*this);
  R2Image tempImage(*this);

  //fprintf(stderr, "applying sobel x and y\n");

  sobelX2.SobelX();
  sobelY2.SobelY();

  //fprintf(stderr, "multiplying x and y into xy\n");

  for (int i=0; i<width; i++) {
    for (int j=0; j<height; j++) {
      sobelXY.Pixel(i,j) = sobelY2.Pixel(i,j) * sobelX2.Pixel(i,j);
    }
  }

  for (int i=0; i<width; i++) {
    for (int j=0; j<height; j++) {
      sobelX2.Pixel(i,j)*= sobelX2.Pixel(i,j);
      sobelY2.Pixel(i,j)*= sobelY2.Pixel(i,j);
    }
  }

  //fprintf(stderr, "applying blur\n");

  sobelX2.Blur(sigma);
  sobelY2.Blur(sigma);
  sobelXY.Blur(sigma);

  //fprintf(stderr, "creating new pixels\n");

  R2Pixel* x2 = new R2Pixel();
  R2Pixel* y2 = new R2Pixel();
  R2Pixel* xy = new R2Pixel();

  //fprintf(stderr, "doing math for loop\n");

  for (int i=0; i<width; i++) {
    for (int j=0; j<height; j++) {

      x2 = &sobelX2.Pixel(i,j);
      y2 = &sobelY2.Pixel(i,j);
      xy = &sobelXY.Pixel(i,j);

      tempImage.Pixel(i,j) = (*x2 * *y2) - (*xy * *xy) - (0.04 * ((*x2 + *y2)) * (*x2 + *y2));
    }
  }
  
  //fprintf(stderr, "creating featureVec\n");

  std::vector<Feature> featureVec;

  R2Pixel* offset = new R2Pixel(.5,.5,.5,1);

//fprintf(stderr, "width: %d height: %d\n", width, height);

  for (int k = 0; k < width; k++) {
    for (int p = 0; p < height; p++) {
        //fprintf(stderr, "starting loop\n");
        tempImage.Pixel(k,p)+= *offset;

        Feature currFeature = Feature(k, p, tempImage.Pixel(k, p));
        featureVec.push_back(currFeature);
        //fprintf(stderr, "pushed back one feature k: %d p: %d \n", k, p);

//***********UNCOMMENT THIS TO HAVE HARRIS FILTER PUT INTO OUTPUT IMAGE***********//

        //Pixel(k,p) = tempImage.Pixel(k,p);

       // fprintf(stderr, "put into output image\n");
    }
  }

  std::sort (featureVec.rbegin(), featureVec.rend());

  // create a double array of variable size
  int variable = this->width*this->height;

  //fprintf(stderr, "variable: %d\n", variable);

  bool * farEnough = (bool*)malloc(variable*sizeof(bool));
  //topFeatureVec.clear();

  for (int i=0; i<height*width; i++) {

    farEnough[i] = false;
  }

  int featuresToMark = 200;
  int featuresMarked = 0;
  int index = 0;
  int buffer = 50;


  while (featuresMarked < featuresToMark) {

    //fprintf(stderr, "while loop\n");

    int currX = featureVec[index].centerX;
    int currY = featureVec[index].centerY;

    //fprintf(stderr, "currX:  %d   currY:  %d \n", currX,currY);

    int indexIn1D = currX + (currY * width); //mapping to 1D array

   // fprintf(stderr, "indexIn1D: %d vs arrayLen: %d\n", indexIn1D, variable);

    bool farEnoughValue = farEnough[indexIn1D];

   // fprintf(stderr, "xvalue: %d\n", farEnoughValue );

    if (farEnoughValue == false) {

      R2Pixel myPixel = Pixel(currX, currY);
      Feature myF(currX, currY, myPixel);
      topFeatureVec.push_back(myF);

      /*for (int j=-(width/256); j<=(width/256); j++) {
        for (int t=-5; t<5; t++){

          Pixel(currX+j+t,currY).SetRed(1);
          Pixel(currX+j+t,currY).SetGreen(0);
          Pixel(currX+j+t,currY).SetBlue(0);

          Pixel(currX,currY+j+t).SetRed(1);
          Pixel(currX,currY+j+t).SetGreen(0);
          Pixel(currX,currY+j+t).SetBlue(0);
      }
    }*/
              
      for (int p= -buffer; p <= buffer; p++) {
        for (int q= -buffer; q <= buffer;q++) {

          int x = currX + p;
          int y = currY + q;
          bool xy_withinRange = (x >= 0 && x < width && y >= 0 && y < height);

          if (xy_withinRange) {
            int temp = (x) + (y * width);

            //fprintf(stderr, "temp: %d vs %d\n", temp, height*width);

            bool temp_withinRange = (temp>=0) && (temp < height*width);

            if (temp_withinRange){
              farEnough[temp] = true;

              // if (abs(p) >= buffer-1 || abs(q) >= buffer-1) {
              // fprintf(stderr, "in if statement:: x: %d y: %d vs width: %d height: %d\n", x, y, width, height);

              // Pixel(x, y).SetRed(1);
              // Pixel(x, y).SetGreen(0);
              // Pixel(x, y).SetBlue(0);
              // fprintf(stderr, "successfully marked box\n");
              // }
            }
          }
        }
      }

      featuresMarked++;
    }
    index++;
  }

  

  // you are supposed to free up pointers at the end
  //printf("trying to end\n");
  free(farEnough);

}


void R2Image::
Sharpen()
{
  // Sharpen an image using a linear filter. Use a kernel of your choosing.
  R2Image* tempImage = new R2Image(width, height);

  R2Pixel outputPixel;

  double kernel[9] = {0.0,-1.0,-0.0,-1.0,5.0,-1.0,0.0,-1.0,0.0}; //kernel is written L>R top>bottom

  for (int i = 1; i < width-1; i++) {
    for (int j = 1; j < height-1; j++) {

      outputPixel.Reset(0.0,0.0,0.0,1.0);

      outputPixel += Pixel(i-1, j+1) * kernel[0];
      outputPixel += Pixel(i, j+1) * kernel[1];
      outputPixel += Pixel(i+1, j+1) * kernel[2];
      outputPixel += Pixel(i-1, j) * kernel[3];
      outputPixel += Pixel(i, j) * kernel[4];
      outputPixel += Pixel(i+1,j) * kernel[5];
      outputPixel += Pixel(i-1,j-1) * kernel[6];
      outputPixel += Pixel(i,j-1) * kernel[7];
      outputPixel += Pixel(i+1,j-1) * kernel[8];

      outputPixel = outputPixel/(kernel[0] + kernel[1] + kernel[2] + kernel[3] + kernel[4] + kernel[5] + kernel[6] + kernel[7] + kernel[8]);

      tempImage->Pixel(i,j) = outputPixel;
      tempImage->Pixel(i,j).Clamp();

    }
  }

  for (int k = 1; k < width-1; k++) {
    for (int p = 1; p < height-1; p++) {

      Pixel(k,p) = tempImage->Pixel(k,p);
    }
  }
}

std::vector<std::vector<double> > markedFeatures;

double** R2Image::
blendOtherImageTranslated(R2Image * otherImage)
{
  // find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
  // compute the matching translation (pixel precision is OK), and blend the translated "otherImage" 
  // into this image with a 50% opacity.
 // fprintf(stderr, "copying *this into tempImage\n");

  R2Image tempImage(*this);

  //fprintf(stderr, "creating currentPixel object\n");

  R2Pixel* currentPixelSSD = new R2Pixel(0,0,0,1);

  tempImage.Harris(2);

  int featuresMarked = 250;
  int boxSideLen = 10;

  for (int i=0; i<featuresMarked; i++) {

    double lowestSum = std::numeric_limits<double>::max();
    int lowestX = 0;
    int lowestY = 0;

    int currX = topFeatureVec[i].centerX;
    int currY = topFeatureVec[i].centerY;

    double currentSum;

    //fprintf(stderr,"accessed topFeatureVec: currX: %d currY: %d\n", currX, currY);

    for (int p=(-width/7+(boxSideLen/2)); p<width/7-(boxSideLen/2); p++) { //p defines window regardless of relative location

      for (int q=(-height/7+(boxSideLen/2)); q<height/7-(boxSideLen/2); q++) {
        //fprintf(stderr, "within window double for loops\n");

        bool withinRange = p+currX >= 0 && p+currX < width && q+currY >= 0 && q+currY < height;

        if (withinRange) {

          currentSum = 0;
          for (int k = -boxSideLen/2; k < boxSideLen/2+1; k++) { //k defines box regardless of relative location

            for (int j= -boxSideLen/2; j<boxSideLen/2+1; j++) {
              //fprintf(stderr, "within box double for loops\n");

              bool withinRangeAgain = (currX + p + k >= 0) && (currX + p + k < std::min(otherImage->width, width)) && (currY + q + j >= 0) && (currY + q + j < std::min(height, otherImage->height));

              if (withinRangeAgain) {

                //fprintf(stderr, "withinRangeAgain\n");

                int u = currX + p + k;
                int v = currY + q + j;

                *currentPixelSSD = (Pixel(currX+k,currY+j)-(otherImage->Pixel(u,v)));

                currentSum += pow(currentPixelSSD->Red(),2);
                currentSum += pow(currentPixelSSD->Green(),2);
                currentSum += pow(currentPixelSSD->Blue(),2);

                //fprintf(stderr, "computed SSD calculations\n");

              }
            }
          }
        }
          if (currentSum < lowestSum) {
            //fprintf(stderr, "currentSum: %f, lowestSum: %f\n", currentSum, lowestSum);

            //fprintf(stderr, "replacing lowestSum lowestX: %d lowestY: %d\n", lowestX, lowestY);

            lowestSum = currentSum;
            lowestX = currX+p;
            lowestY = currY+q;
          }
        }
      }
      

      //fprintf(stderr, "currX: %d currY %d lowestX: %d, lowestY: %d", currX, currY, lowestX, lowestY);
      
    topFeatureVec[i].endX = lowestX;
    topFeatureVec[i].endY = lowestY;

      //**********tempImage2.line(currX,lowestX,currY,lowestY, 1.0, 0.0, 0.0);

    }
int iterations = 1000;
int maxMatches = 0;
int currentMatches = 0;
int threshhold = 7;
double** bestHMatrix = dmatrix(1,3,1,3);
double** perfHMatrix = dmatrix(1,3,1,3);
std::vector<double> goodFeatures;
srand(time(NULL));

for (int iteration = 0; iteration<iterations; iteration++) {

  currentMatches = 0;

    int f1Index = rand()%(featuresMarked);
    int f2Index = rand()%(featuresMarked);
    int f3Index = rand()%(featuresMarked);
    int f4Index = rand()%(featuresMarked);

    Feature tempFeature = topFeatureVec[f1Index];

    double x1 = tempFeature.centerX;
    double y1 = tempFeature.centerY;
    double x2 = tempFeature.endX;
    double y2 = tempFeature.endY;

    fprintf(stderr, "FEATURE 1: Index = %d \n x1 = %f y1 = %f \n x2 = %f y2 = %f \n",f1Index, x1, y1, x2, y2);

    tempFeature = topFeatureVec[f2Index];

    double x3 = tempFeature.centerX;
    double y3 = tempFeature.centerY;
    double x4 = tempFeature.endX;
    double y4 = tempFeature.endY;

    fprintf(stderr, "FEATURE 2: Index = %d \n x3 = %f y3 = %f \n x4 = %f y4 = %f \n",f2Index, x3, y3, x4, y4);


    tempFeature = topFeatureVec[f3Index];

    double x5 = tempFeature.centerX;
    double y5 = tempFeature.centerY;
    double x6 = tempFeature.endX;
    double y6 = tempFeature.endY;

    fprintf(stderr, "FEATURE 3: Index = %d \n x5 = %f y5 = %f \n x6 = %f y6 = %f \n",f3Index, x5, y5, x6, y6);

    tempFeature = topFeatureVec[f4Index];

    double x7 = tempFeature.centerX;
    double y7 = tempFeature.centerY;
    double x8 = tempFeature.endX;
    double y8 = tempFeature.endY;

    fprintf(stderr, "FEATURE 4: Index = %d \n x7 = %f y7 = %f \n x8 = %f y8 = %f \n",f4Index, x7, y7, x8, y8);

    std::vector<double> input;

    input.push_back(x1);
    input.push_back(y1);
    input.push_back(x2);
    input.push_back(y2);
    input.push_back(x3);
    input.push_back(y3);
    input.push_back(x4);
    input.push_back(y4);
    input.push_back(x5);
    input.push_back(y5);
    input.push_back(x6);
    input.push_back(y6);
    input.push_back(x7);
    input.push_back(y7);
    input.push_back(x8);
    input.push_back(y8);

    //fprintf(stderr, "input Size: %lu\n", input.size());
    double** hMatrix = svdTest(input);

    //printf("CURRENT H matrix: \n%f %f %f\n %f %f %f\n %f %f %f\n", hMatrix[1][1], hMatrix[1][2], hMatrix[1][3], hMatrix[2][1], hMatrix[2][2], hMatrix[2][3], hMatrix[3][1], hMatrix[3][2], hMatrix[3][3]);


    for (int i=0; i<featuresMarked; i++) {
      if (i!=f1Index && i!=f2Index && i!=f3Index && i!=f4Index) {

        double originVecEndX = topFeatureVec[i].endX;
        double originVecEndY = topFeatureVec[i].endY;
        double originVecCenterX = topFeatureVec[i].centerX;
        double originVecCenterY = topFeatureVec[i].centerY;

        double expectedVecEndZ = ((originVecCenterX)*(hMatrix[3][1])) + ((originVecCenterY)*(hMatrix[3][2])) + (hMatrix[3][3]);

        double expectedVecEndX = ((originVecCenterX)*(hMatrix[1][1]) + (originVecCenterY)*(hMatrix[1][2]) + (hMatrix[1][3]))/expectedVecEndZ;
        double expectedVecEndY = ((originVecCenterX)*(hMatrix[2][1]) + (originVecCenterY)*(hMatrix[2][2]) + (hMatrix[2][3]))/expectedVecEndZ;

        fprintf(stderr, "expectedVecEndZ = (%f)*(%f) + (%f)*(%f) + (%f) = %f\n", originVecCenterX, bestHMatrix[3][1],originVecCenterY, bestHMatrix[3][2], bestHMatrix[3][3], expectedVecEndZ);
        fprintf(stderr, "expectedVecEndX = ((%f)*(%f) + (%f)*(%f) + (%f))/%f = %f\n", originVecCenterX, bestHMatrix[3][1],originVecCenterY, bestHMatrix[3][2], bestHMatrix[3][3], expectedVecEndZ,expectedVecEndX);
        fprintf(stderr, "expectedVecEndY = ((%f)*(%f) + (%f)*(%f) + (%f))/%f = %f\n", originVecCenterX, bestHMatrix[3][1],originVecCenterY, bestHMatrix[3][2], bestHMatrix[3][3], expectedVecEndZ,expectedVecEndY);

        fprintf(stderr, "distance(%f, %f, %f, %f)\n",originVecEndX, expectedVecEndX, originVecEndY, expectedVecEndY);
        double temp = distance(originVecEndX, expectedVecEndX, originVecEndY, expectedVecEndY);

        fprintf(stderr, "i = %d; temp = %f\n",i, temp);

        if (temp <= threshhold) {
          //fprintf(stderr, "if\n" );
          currentMatches ++;
        }
      }
    }

fprintf(stderr, "currentMatches = %d\n", currentMatches);
    if (currentMatches > maxMatches) {

      maxMatches = currentMatches;

      for (int q=1; q<4; q++) {
        for (int p=1; p<4; p++) {
          bestHMatrix[p][q] = hMatrix[p][q];
        }
      }
    }
  }

  R2Image tempImage2(*this);

  printf("BEST H matrix: \n%f %f %f\n %f %f %f\n %f %f %f\n", bestHMatrix[1][1], bestHMatrix[1][2], bestHMatrix[1][3], bestHMatrix[2][1], bestHMatrix[2][2], bestHMatrix[2][3], bestHMatrix[3][1], bestHMatrix[3][2], bestHMatrix[3][3]);

  for (int i=0; i<250; i++) {

    double originVecEndX = topFeatureVec[i].endX;
    double originVecEndY = topFeatureVec[i].endY;
    double originVecCenterX = topFeatureVec[i].centerX;
    double originVecCenterY = topFeatureVec[i].centerY;

    // markedFeatures[i][0] = originVecCenterX;
    // markedFeatures[i][1] = originVecCenterY;
    // markedFeatures[i][2] = originVecEndX;
    // markedFeatures[i][3] = originVecEndY;

    double expectedVecEndZ = ((originVecCenterX)*(bestHMatrix[3][1])) + ((originVecCenterY)*(bestHMatrix[3][2])) + (bestHMatrix[3][3]);

    double expectedVecEndX = (((originVecCenterX)*(bestHMatrix[1][1])) + ((originVecCenterY)*(bestHMatrix[1][2])) + (bestHMatrix[1][3]))/expectedVecEndZ;
    double expectedVecEndY = (((originVecCenterX)*(bestHMatrix[2][1])) + ((originVecCenterY)*(bestHMatrix[2][2])) + (bestHMatrix[2][3]))/expectedVecEndZ;

    double temp = distance(originVecEndX, expectedVecEndX, originVecEndY, expectedVecEndY);


    //fprintf(stderr, "i = %d; temp = %f\n",i, temp);

    if (temp <= threshhold) {
      goodFeatures.push_back(originVecCenterX);
      goodFeatures.push_back(originVecCenterY);
      goodFeatures.push_back(originVecEndX);
      goodFeatures.push_back(originVecEndY);

      // markedFeatures[i][4] = 1;

      tempImage2.line(originVecCenterX,originVecEndX,originVecCenterY,originVecEndY, 0.0, 1.0, 0.0);

      tempImage2.line(originVecCenterX, expectedVecEndX, originVecCenterY, expectedVecEndY, 0.541176,0.168627, 0.886275);

    }
    else {
      // markedFeatures[i][4] = 0;
      tempImage2.line(originVecCenterX,originVecEndX,originVecCenterY,originVecEndY, 1.0, 0.0, 0.0);
    }
  }

  perfHMatrix = svdTest(goodFeatures);
  printf("PERF H matrix: \n%f %f %f\n %f %f %f\n %f %f %f\n", perfHMatrix[1][1], perfHMatrix[1][2], perfHMatrix[1][3], perfHMatrix[2][1], perfHMatrix[2][2], perfHMatrix[2][3], perfHMatrix[3][1], perfHMatrix[3][2], perfHMatrix[3][3]);

fprintf(stderr, "1011\n");
    for (int p=0; p<width; p++) {
      for (int q=0; q<height; q++) {
        Pixel(p,q) = tempImage2.Pixel(p,q);
      }
    }
    fprintf(stderr, "1017\n");
  return perfHMatrix;
  }


double R2Image::distance(int x1, int x2, int y1, int y2) {

  double result = sqrt(pow((x2-x1),2) + pow((y2-y1),2));
  return result;
}

void R2Image::
blendOtherImageHomography(R2Image * otherImage)
{

  R2Image warpedSecondImage(*this);
  R2Image copyFirstImage(*this);
  R2Pixel pixelToCopy = R2black_pixel;
  R2Pixel pixelToCopy2 = R2black_pixel;
  double** bestHMatrix = blendOtherImageTranslated(otherImage);

  // fprintf(stderr, "%s\n", "1245");
  for (int i=0; i<width; i++) {

    for (int j=0; j<height; j++) {
      double expectedVecEndZ = ((i)*(bestHMatrix[3][1])) + ((j)*(bestHMatrix[3][2])) + (bestHMatrix[3][3]);

      double expectedVecEndX = (((i)*(bestHMatrix[1][1])) + ((j)*(bestHMatrix[1][2])) + (bestHMatrix[1][3]))/expectedVecEndZ;
      double expectedVecEndY = (((i)*(bestHMatrix[2][1])) + ((j)*(bestHMatrix[2][2])) + (bestHMatrix[2][3]))/expectedVecEndZ;

      // fprintf(stderr, "projZ: %f projX: %f projY: %f\n", expectedVecEndZ, expectedVecEndX, expectedVecEndY);
      pixelToCopy = otherImage->Pixel(expectedVecEndX, expectedVecEndY);
      warpedSecondImage.Pixel(i,j) = pixelToCopy;
    }
  }

  for (int p=0; p<width; p++) {
    for (int q=0; q<height; q++) {

      pixelToCopy2 = warpedSecondImage.Pixel(p,q) + this->Pixel(p,q);
      pixelToCopy2 /= 2;

      copyFirstImage.Pixel(p,q) = pixelToCopy2;
    }
  }

  for (int j=0; j<width; j++) {
    for (int k=0; k<height; k++) {

      Pixel(j,k) = copyFirstImage.Pixel(j,k);
    }
  }
}

R2Image* copyImage;
double averageDeltaX;
double averageDeltaY;

void R2Image::
  firstFrameProcessing() {

    copyImage = new R2Image (*this);

    this->Harris(3);

    averageDeltaY = 0;
    averageDeltaX = 0;

    //Create copy of topFeatureVec as featuresToCheck
    int featuresMarked = 200;
    numFeaturesToCheck = 200;
    for(int i = 0; i < featuresMarked; i++){
      featuresToCheck.push_back(topFeatureVec[i]);
    }

    return;
}



void R2Image::
  frameProcessing(R2Image* otherImage) {

    //numFeaturesToCheck--;
    fprintf(stderr, "%d\n", numFeaturesToCheck);
    //numFeaturesToCheck= numFeaturesToCheck -1;
    //fprintf(stderr, "%d\n", numFeaturesToCheck);


  R2Image tempImage(*this);

  // fprintf(stderr, "creating currentPixel object\n");

  R2Pixel* currentPixelSSD = new R2Pixel(0,0,0,1);

  //int featuresMarked = 250;
  int boxSideLen = 7;

// fprintf(stderr, "frameProcessing\n" );

  for (int i=0; i<numFeaturesToCheck; i++) {

    double lowestSum = std::numeric_limits<double>::max();
    int lowestX = 0;
    int lowestY = 0;

    int currX = featuresToCheck[i].centerX;
    int currY = featuresToCheck[i].centerY;

    double currentSum;

    //fprintf(stderr,"accessed topFeatureVec: currX: %d currY: %d\n", currX, currY);

    for (int p=(-width/35); p<width/35; p++) { //p defines window regardless of relative location

      for (int q=(-height/35); q<height/35; q++) {
     
    // for (int p=(-width/25+(boxSideLen/2)); p<width/25-(boxSideLen/2); p++) { //p defines window regardless of relative location

    //   for (int q=(-height/25+(boxSideLen/2)); q<height/25-(boxSideLen/2); q++) {
    //fprintf(stderr, "within window double for loops\n");

        bool withinRange = p+currX >= 0 && p+currX < width && q+currY >= 0 && q+currY < height;

        if (withinRange) {

          currentSum = 0;
          for (int k = -boxSideLen/2; k < boxSideLen/2+1; k++) { //k defines box regardless of relative location

            for (int j= -boxSideLen/2; j<boxSideLen/2+1; j++) {
              //fprintf(stderr, "within box double for loops\n");

              bool withinRangeAgain = (currX + p + k >= 0) && (currX + p + k < std::min(otherImage->width, width)) && (currY + q + j >= 0) && (currY + q + j < std::min(height, otherImage->height));

              if (withinRangeAgain) {

                //fprintf(stderr, "withinRangeAgain\n");

                int u = currX + p + k;
                int v = currY + q + j;

                *currentPixelSSD = (Pixel(currX+k,currY+j)-(otherImage->Pixel(u,v)));


                currentSum += pow(currentPixelSSD->Red(),2);
                currentSum += pow(currentPixelSSD->Green(),2);
                currentSum += pow(currentPixelSSD->Blue(),2);

              }
            }
          }
        }
          if (currentSum < lowestSum) {
            //fprintf(stderr, "currentSum: %f, lowestSum: %f\n", currentSum, lowestSum);

            //fprintf(stderr, "replacing lowestSum lowestX: %d lowestY: %d\n", lowestX, lowestY);

            lowestSum = currentSum;
            lowestX = currX+p;
            lowestY = currY+q;
          }
      }
    }
      

      //fprintf(stderr, "currX: %d currY %d lowestX: %d, lowestY: %d", currX, currY, lowestX, lowestY);
      
    featuresToCheck[i].endX = lowestX;
    featuresToCheck[i].endY = lowestY;

      //**********tempImage2.line(currX,lowestX,currY,lowestY, 1.0, 0.0, 0.0);

  }

  //BEGIN SIMPLIFIED RANSAC WITH 1 RANDOM FEATURE (10 ITERATIONS)
  int maxSize = 0;
  int maxIndex;
  int currentSize;
  double threshhold = 8;
  for (int i=0; i<8; i++) {

    currentSize = 0;
    int randomNum = rand()%numFeaturesToCheck;
    int matchVecX1 = featuresToCheck[randomNum].centerX; // 5
    int matchVecY1 = featuresToCheck[randomNum].centerY; // 1
    int matchVecX2 = featuresToCheck[randomNum].endX; // 8
    int matchVecY2 = featuresToCheck[randomNum].endY; // 2

    for (int j=0; j<numFeaturesToCheck; j++) {

        int compareVecX1 = featuresToCheck[j].centerX; // 1
        int compareVecY1 = featuresToCheck[j].centerY; // 1
        int compareVecX2 = featuresToCheck[j].endX; // 3
        int compareVecY2 = featuresToCheck[j].endY; // 2

        int offsetX = compareVecX1-matchVecX1; //1-5 = -4
        int offsetY = compareVecY1-matchVecY1; // 1-1 = 0


        double temp = distance(compareVecX2, matchVecX2+offsetX, compareVecY2, matchVecY2+offsetY);

        if (temp <= threshhold) {

          currentSize ++;
        }
      }
    if (currentSize > maxSize) {
        maxSize = currentSize;
        maxIndex = randomNum;
      }    
  }

  //DO RANSAC AGAIN WITH BEST MATCH VECTOR

  int matchVecX1 = featuresToCheck[maxIndex].centerX; // 5
  int matchVecY1 = featuresToCheck[maxIndex].centerY; // 1
  int matchVecX2 = featuresToCheck[maxIndex].endX; // 8
  int matchVecY2 = featuresToCheck[maxIndex].endY; // 2

  std::vector<double> deleteIndexList;
  double avgDeltaXFrame = 0;
  double avgDeltaYFrame = 0;

  int inliers = 0;
  for (int i=0; i<numFeaturesToCheck; i++) {
    // tempFeatureList.clear();
    //deleteIndexList.clear();
    int compareVecX1 = featuresToCheck[i].centerX; // 1
    int compareVecY1 = featuresToCheck[i].centerY; // 1
    int compareVecX2 = featuresToCheck[i].endX; // 3
    int compareVecY2 = featuresToCheck[i].endY; // 2

    int offsetX = compareVecX1-matchVecX1; //1-5 = -4
    int offsetY = compareVecY1-matchVecY1; // 1-1 = 0

    double temp = distance(compareVecX2, matchVecX2+offsetX, compareVecY2, matchVecY2+offsetY);

    if (temp <= threshhold) {
      tempImage.line(compareVecX1,compareVecX2,compareVecY1,compareVecY2, 0.0, 1.0, 0.0);
      avgDeltaXFrame += double(compareVecX2 - compareVecX1);
      avgDeltaYFrame += double(compareVecY2 - compareVecY1);
      inliers++;
    }
    else {
      tempImage.line(compareVecX1,compareVecX2,compareVecY1,compareVecY2, 1.0, 0.0, 0.0);
      //Add bad feature index to delete index list
      deleteIndexList.push_back(i);
    }
  //Reset feature centers as ends for next frame pair after marking feature
  featuresToCheck[i].centerX = featuresToCheck[i].endX;
  featuresToCheck[i].centerY = featuresToCheck[i].endY;
  }

  avgDeltaXFrame = avgDeltaXFrame / inliers;
  avgDeltaYFrame = avgDeltaYFrame / inliers;

  //printf("PERF H matrix: \n%f %f %f\n %f %f %f\n %f %f %f\n", perfHMatrix[1][1], perfHMatrix[1][2], perfHMatrix[1][3], perfHMatrix[2][1], perfHMatrix[2][2], perfHMatrix[2][3], perfHMatrix[3][1], perfHMatrix[3][2], perfHMatrix[3][3]);
      //std::vector<double> goodFeatures;
    
  //Delete bad features from featuresToCheck list
  int idx;
  int deleteIndexListSize = static_cast<int>(deleteIndexList.size());
  //fprintf(stderr, "size of featuresToCheck: %d\n", numFeaturesToCheck);
  //fprintf(stderr, "size of deleteIndexList: %d\n" ,deleteIndexListSize-1);

  for(int i=deleteIndexListSize-1; i > 0; i--) {
      
      idx = deleteIndexList[i];

      numFeaturesToCheck = numFeaturesToCheck -1;
      
      for (int j=idx; j<numFeaturesToCheck; j++) {
        featuresToCheck[j] = featuresToCheck[j+1];

      }
  }

  deleteIndexList.clear();

  averageDeltaX += avgDeltaXFrame;
  averageDeltaY += avgDeltaYFrame;


    for (int p=0; p<width; p++) {
      for (int q=0; q<height; q++) {
        Pixel(p,q) = tempImage.Pixel(p,q);
      }
    }
}

double** R2Image::computeStabilizationMatrix() {

  averageDeltaX /= 88;
  averageDeltaY /= 88;

  double** result = dmatrix(1,3,1,3);

  result[1][1] = 1;
  result[1][2] = 0;
  result[1][3] = averageDeltaX;
  result[2][1] = 0;
  result[2][2] = 1;
  result[2][3] = averageDeltaY;
  result[3][1] = 0;
  result[3][2] = 0;
  result[3][3] = 1;

  return result;
}

void R2Image::stabilization(R2Image* otherImage, double** stabilizationMatrix) {


  R2Image *copyImage = new R2Image(*this);
    fprintf(stderr, "copy image created\n" );


  for (int i=0; i< width; i++) {
    for (int j=0; j<height; j++) {
      double expectedEndZ = ((i)*(stabilizationMatrix[3][1])) + ((j)*(stabilizationMatrix[3][2])) + (stabilizationMatrix[3][3]);

      double expectedEndX = (((i)*(stabilizationMatrix[1][1])) + ((j)*(stabilizationMatrix[1][2])) + (stabilizationMatrix[1][3]))/expectedEndZ;
      double expectedEndY = (((i)*(stabilizationMatrix[2][1])) + ((j)*(stabilizationMatrix[2][2])) + (stabilizationMatrix[2][3]))/expectedEndZ;

      fprintf(stderr, "copyImage->Pixel(%f, %f) = this->Pixel(%d,%d);\n",expectedEndX, expectedEndY, i, j );

      bool validPoint = expectedEndX > 0 && expectedEndY > 0 && expectedEndX < (width-1) && expectedEndY < (height-1);
      if (validPoint) {
      copyImage->Pixel(expectedEndX, expectedEndY) = this->Pixel(i,j);
      }
    }
  }

  for (int i=0; i< width; i++) {
    for (int j=0; j<height; j++) {

      Pixel(i,j) = copyImage->Pixel(i,j);
    }
  }

}

// void R2Image::
//   frameProcessing(R2Image* otherImage) {

//     //numFeaturesToCheck--;
//     fprintf(stderr, "%d\n", numFeaturesToCheck);
//     //numFeaturesToCheck= numFeaturesToCheck -1;
//     //fprintf(stderr, "%d\n", numFeaturesToCheck);


//   R2Image tempImage(*this);

//   // fprintf(stderr, "creating currentPixel object\n");

//   R2Pixel* currentPixelSSD = new R2Pixel(0,0,0,1);

//   //int featuresMarked = 250;
//   int boxSideLen = 8;

// // fprintf(stderr, "frameProcessing\n" );

//   for (int i=0; i<numFeaturesToCheck; i++) {

//     double lowestSum = std::numeric_limits<double>::max();
//     int lowestX = 0;
//     int lowestY = 0;

//     int currX = featuresToCheck[i].centerX;
//     int currY = featuresToCheck[i].centerY;

//     double currentSum;

//     //fprintf(stderr,"accessed topFeatureVec: currX: %d currY: %d\n", currX, currY);

// for (int p=(-width/25); p<width/25; p++) { //p defines window regardless of relative location

//       for (int q=(-height/25); q<height/25; q++) {
     
//     // for (int p=(-width/25+(boxSideLen/2)); p<width/25-(boxSideLen/2); p++) { //p defines window regardless of relative location

//     //   for (int q=(-height/25+(boxSideLen/2)); q<height/25-(boxSideLen/2); q++) {
//     //fprintf(stderr, "within window double for loops\n");

//         bool withinRange = p+currX >= 0 && p+currX < width && q+currY >= 0 && q+currY < height;

//         if (withinRange) {

//           currentSum = 0;
//           for (int k = -boxSideLen/2; k < boxSideLen/2+1; k++) { //k defines box regardless of relative location

//             for (int j= -boxSideLen/2; j<boxSideLen/2+1; j++) {
//               //fprintf(stderr, "within box double for loops\n");

//               bool withinRangeAgain = (currX + p + k >= 0) && (currX + p + k < std::min(otherImage->width, width)) && (currY + q + j >= 0) && (currY + q + j < std::min(height, otherImage->height));

//               if (withinRangeAgain) {

//                 //fprintf(stderr, "withinRangeAgain\n");

//                 int u = currX + p + k;
//                 int v = currY + q + j;

//                 *currentPixelSSD = (Pixel(currX+k,currY+j)-(otherImage->Pixel(u,v)));


//                 currentSum += pow(currentPixelSSD->Red(),2);
//                 currentSum += pow(currentPixelSSD->Green(),2);
//                 currentSum += pow(currentPixelSSD->Blue(),2);

//               }
//             }
//           }
//         }
//           if (currentSum < lowestSum) {
//             //fprintf(stderr, "currentSum: %f, lowestSum: %f\n", currentSum, lowestSum);

//             //fprintf(stderr, "replacing lowestSum lowestX: %d lowestY: %d\n", lowestX, lowestY);

//             lowestSum = currentSum;
//             lowestX = currX+p;
//             lowestY = currY+q;
//           }
//         }
//       }
      

//       //fprintf(stderr, "currX: %d currY %d lowestX: %d, lowestY: %d", currX, currY, lowestX, lowestY);
      
//     featuresToCheck[i].endX = lowestX;
//     featuresToCheck[i].endY = lowestY;

//       //**********tempImage2.line(currX,lowestX,currY,lowestY, 1.0, 0.0, 0.0);

//     }
// int iterations = 1;
// int maxMatches = 0;
// int currentMatches = 0;
// int threshhold = 20;
// double** bestHMatrix = dmatrix(1,3,1,3);
// double** perfHMatrix = dmatrix(1,3,1,3);
// std::vector<double> goodFeatures;
// srand(time(NULL));

// for (int iteration = 0; iteration<iterations; iteration++) {

//     currentMatches = 0;

//     int f1Index = rand()%(numFeaturesToCheck);
//     int f2Index = rand()%(numFeaturesToCheck);
//     int f3Index = rand()%(numFeaturesToCheck);
//     int f4Index = rand()%(numFeaturesToCheck);

//     Feature tempFeature = featuresToCheck[f1Index];

//     double x1 = tempFeature.centerX;
//     double y1 = tempFeature.centerY;
//     double x2 = tempFeature.endX;
//     double y2 = tempFeature.endY;

//     fprintf(stderr, "FEATURE 1: Index = %d \n x1 = %f y1 = %f \n x2 = %f y2 = %f \n",f1Index, x1, y1, x2, y2);

//     tempFeature = featuresToCheck[f2Index];

//     double x3 = tempFeature.centerX;
//     double y3 = tempFeature.centerY;
//     double x4 = tempFeature.endX;
//     double y4 = tempFeature.endY;

//     fprintf(stderr, "FEATURE 2: Index = %d \n x3 = %f y3 = %f \n x4 = %f y4 = %f \n",f2Index, x3, y3, x4, y4);

//     tempFeature = featuresToCheck[f3Index];

//     double x5 = tempFeature.centerX;
//     double y5 = tempFeature.centerY;
//     double x6 = tempFeature.endX;
//     double y6 = tempFeature.endY;

//     fprintf(stderr, "FEATURE 3: Index = %d \n x5 = %f y5 = %f \n x6 = %f y6 = %f \n",f3Index, x5, y5, x6, y6);

//     tempFeature = featuresToCheck[f4Index];

//     double x7 = tempFeature.centerX;
//     double y7 = tempFeature.centerY;
//     double x8 = tempFeature.endX;
//     double y8 = tempFeature.endY;

//     fprintf(stderr, "FEATURE 4: Index = %d \n x7 = %f y7 = %f \n x8 = %f y8 = %f \n",f4Index, x7, y7, x8, y8);

//     std::vector<double> input;

//     input.push_back(x1);
//     input.push_back(y1);
//     input.push_back(x2);
//     input.push_back(y2);
//     input.push_back(x3);
//     input.push_back(y3);
//     input.push_back(x4);
//     input.push_back(y4);
//     input.push_back(x5);
//     input.push_back(y5);
//     input.push_back(x6);
//     input.push_back(y6);
//     input.push_back(x7);
//     input.push_back(y7);
//     input.push_back(x8);
//     input.push_back(y8);

//     //fprintf(stderr, "input Size: %lu\n", input.size());
//     double** hMatrix = svdTest(input);

//     // printf("CURRENT H matrix: \n%f %f %f\n %f %f %f\n %f %f %f\n", hMatrix[1][1], hMatrix[1][2], hMatrix[1][3], hMatrix[2][1], hMatrix[2][2], hMatrix[2][3], hMatrix[3][1], hMatrix[3][2], hMatrix[3][3]);

//     for (int i=0; i<numFeaturesToCheck; i++) {
//       //if (i!=f1Index && i!=f2Index && i!=f3Index && i!=f4Index) {

//         double originVecEndX = featuresToCheck[i].endX;
//         double originVecEndY = featuresToCheck[i].endY;
//         double originVecCenterX = featuresToCheck[i].centerX;
//         double originVecCenterY = featuresToCheck[i].centerY;

//         double expectedVecEndZ = ((originVecCenterX)*(hMatrix[3][1])) + ((originVecCenterY)*(hMatrix[3][2])) + (hMatrix[3][3]);

//         double expectedVecEndX = ((originVecCenterX)*(hMatrix[1][1]) + (originVecCenterY)*(hMatrix[1][2]) + (hMatrix[1][3]))/expectedVecEndZ;
//         double expectedVecEndY = ((originVecCenterX)*(hMatrix[2][1]) + (originVecCenterY)*(hMatrix[2][2]) + (hMatrix[2][3]))/expectedVecEndZ;

//         // fprintf(stderr, "expectedVecEndZ = (%f)*(%f) + (%f)*(%f) + (%f) = %f\n", originVecCenterX, bestHMatrix[3][1],originVecCenterY, bestHMatrix[3][2], bestHMatrix[3][3], expectedVecEndZ);
//         // fprintf(stderr, "expectedVecEndX = ((%f)*(%f) + (%f)*(%f) + (%f))/%f = %f\n", originVecCenterX, bestHMatrix[3][1],originVecCenterY, bestHMatrix[3][2], bestHMatrix[3][3], expectedVecEndZ,expectedVecEndX);
//         // fprintf(stderr, "expectedVecEndY = ((%f)*(%f) + (%f)*(%f) + (%f))/%f = %f\n", originVecCenterX, bestHMatrix[3][1],originVecCenterY, bestHMatrix[3][2], bestHMatrix[3][3], expectedVecEndZ,expectedVecEndY);

//         double temp = distance(originVecEndX, expectedVecEndX, originVecEndY, expectedVecEndY);

//         if (i==f1Index || i == f2Index || i == f3Index || i == f4Index) {
//         fprintf(stderr, "distance(%f, %f, %f, %f) = %f\n",originVecEndX, expectedVecEndX, originVecEndY, expectedVecEndY, temp);
//         }
//         //fprintf(stderr, "i = %d; temp = %f\n",i, temp);

//         if (temp <= threshhold) {
//           fprintf(stderr, "within Threshold\n" );
//           currentMatches ++;
//         }
//      // }
//     }

// //fprintf(stderr, "currentMatches = %d\n", currentMatches);
//     if (currentMatches > maxMatches) {

//       maxMatches = currentMatches;

//       for (int q=1; q<4; q++) {
//         for (int p=1; p<4; p++) {
//           bestHMatrix[p][q] = hMatrix[p][q];
//         }
//       }
//     }
//   }

//   R2Image tempImage2(*this);

//   //printf("BEST H matrix: \n%f %f %f\n %f %f %f\n %f %f %f\n", bestHMatrix[1][1], bestHMatrix[1][2], bestHMatrix[1][3], bestHMatrix[2][1], bestHMatrix[2][2], bestHMatrix[2][3], bestHMatrix[3][1], bestHMatrix[3][2], bestHMatrix[3][3]);

//   for (int i=0; i<250; i++) {

//     double originVecEndX = featuresToCheck[i].endX;
//     double originVecEndY = featuresToCheck[i].endY;
//     double originVecCenterX = featuresToCheck[i].centerX;
//     double originVecCenterY = featuresToCheck[i].centerY;

//     // markedFeatures[i][0] = originVecCenterX;
//     // markedFeatures[i][1] = originVecCenterY;
//     // markedFeatures[i][2] = originVecEndX;
//     // markedFeatures[i][3] = originVecEndY;

//     double expectedVecEndZ = ((originVecCenterX)*(bestHMatrix[3][1])) + ((originVecCenterY)*(bestHMatrix[3][2])) + (bestHMatrix[3][3]);

//     double expectedVecEndX = (((originVecCenterX)*(bestHMatrix[1][1])) + ((originVecCenterY)*(bestHMatrix[1][2])) + (bestHMatrix[1][3]))/expectedVecEndZ;
//     double expectedVecEndY = (((originVecCenterX)*(bestHMatrix[2][1])) + ((originVecCenterY)*(bestHMatrix[2][2])) + (bestHMatrix[2][3]))/expectedVecEndZ;

//     double temp = distance(originVecEndX, expectedVecEndX, originVecEndY, expectedVecEndY);


//     //fprintf(stderr, "i = %d; temp = %f\n",i, temp);

//     if (temp <= threshhold) {
//      // fprintf(stderr, "withinThreshold\n");
//       goodFeatures.push_back(originVecCenterX);
//       goodFeatures.push_back(originVecCenterY);
//       goodFeatures.push_back(originVecEndX);
//       goodFeatures.push_back(originVecEndY);

//       // markedFeatures[i][4] = 1;

//       //tempImage2.line(originVecCenterX,originVecEndX,originVecCenterY,originVecEndY, 0.0, 1.0, 0.0);

//       //tempImage2.line(originVecCenterX, expectedVecEndX, originVecCenterY, expectedVecEndY, 0.541176,0.168627, 0.886275);

//     }
//     else {
//       // markedFeatures[i][4] = 0;
//       //tempImage2.line(originVecCenterX,originVecEndX,originVecCenterY,originVecEndY, 1.0, 0.0, 0.0);
//     }
//   }

//   perfHMatrix = svdTest(goodFeatures);
//   //std::vector<Feature> tempFeatureList;
//   std::vector<double> deleteIndexList;
//   //printf("PERF H matrix: \n%f %f %f\n %f %f %f\n %f %f %f\n", perfHMatrix[1][1], perfHMatrix[1][2], perfHMatrix[1][3], perfHMatrix[2][1], perfHMatrix[2][2], perfHMatrix[2][3], perfHMatrix[3][1], perfHMatrix[3][2], perfHMatrix[3][3]);
//       //std::vector<double> goodFeatures;

//       for (int i=0; i<numFeaturesToCheck; i++) {
//        // tempFeatureList.clear();
//       //deleteIndexList.clear();
//       double originVecEndX = featuresToCheck[i].endX;
//       double originVecEndY = featuresToCheck[i].endY;
//       double originVecCenterX = featuresToCheck[i].centerX;
//       double originVecCenterY = featuresToCheck[i].centerY;
  
//       double expectedVecEndZperf = ((originVecCenterX)*(perfHMatrix[3][1])) + ((originVecCenterY)*(perfHMatrix[3][2])) + (perfHMatrix[3][3]);

//       double expectedVecEndXperf = (((originVecCenterX)*(perfHMatrix[1][1])) + ((originVecCenterY)*(perfHMatrix[1][2])) + (perfHMatrix[1][3]))/expectedVecEndZperf;
//       double expectedVecEndYperf = (((originVecCenterX)*(perfHMatrix[2][1])) + ((originVecCenterY)*(perfHMatrix[2][2])) + (perfHMatrix[2][3]))/expectedVecEndZperf;


//       double temp = distance(originVecEndX, expectedVecEndXperf, originVecEndY, expectedVecEndYperf);

//     //fprintf(stderr, "i = %d; temp = %f\n",i, temp);

//       if (temp <= threshhold) {
//        // fprintf(stderr, "withinThreshold\n");

//         // markedFeatures[i][4] = 1;

//         tempImage2.line(originVecCenterX,originVecEndX,originVecCenterY,originVecEndY, 0.0, 1.0, 0.0);

//         //tempImage2.line(originVecCenterX, expectedVecEndX, originVecCenterY, expectedVecEndY, 0.541176,0.168627, 0.886275);

//       }
//       else {
//        // fprintf(stderr, "else\n" );
//         // markedFeatures[i][4] = 0;
//         tempImage2.line(originVecCenterX,originVecEndX,originVecCenterY,originVecEndY, 1.0, 0.0, 0.0);
//         //Add bad feature index to delete index list
//         deleteIndexList.push_back(i);
//     }

//     featuresToCheck[i].centerX = featuresToCheck[i].endX;
//     featuresToCheck[i].centerY = featuresToCheck[i].endY;
//   }
    
//   //Delete bad features from featuresToCheck list
//   int idx;
//   int deleteIndexListSize = static_cast<int>(deleteIndexList.size());
//   //fprintf(stderr, "size of featuresToCheck: %d\n", numFeaturesToCheck);
//   //fprintf(stderr, "size of deleteIndexList: %d\n" ,deleteIndexListSize-1);

//   for(int i=deleteIndexListSize-1; i > 0; i--) {
      
//       idx = deleteIndexList[i];

//       numFeaturesToCheck = numFeaturesToCheck -1;
      
//       for (int j=idx; j<numFeaturesToCheck; j++) {
//         featuresToCheck[j] = featuresToCheck[j+1];

//       }
//   }

//   deleteIndexList.clear();


//     for (int p=0; p<width; p++) {
//       for (int q=0; q<height; q++) {
//         Pixel(p,q) = tempImage2.Pixel(p,q);
//       }
//     }
// }

////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

int R2Image::
Read(const char *filename)
{
  // Initialize everything
  if (pixels) { delete [] pixels; pixels = NULL; }
  npixels = width = height = 0;

  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }
  
  // Read file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return ReadBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return ReadPPM(filename);
  else if (!strncmp(input_extension, ".jpg", 4)) return ReadJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return ReadJPEG(filename);
  
  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



int R2Image::
Write(const char *filename) const
{
  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }
  
  // Write file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return WriteBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return WritePPM(filename, 1);
  else if (!strncmp(input_extension, ".jpg", 5)) return WriteJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return WriteJPEG(filename);

  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



////////////////////////////////////////////////////////////////////////
// BMP I/O
////////////////////////////////////////////////////////////////////////

#if (RN_OS == RN_LINUX) && !WIN32

typedef struct tagBITMAPFILEHEADER {
  unsigned short int bfType;
  unsigned int bfSize;
  unsigned short int bfReserved1;
  unsigned short int bfReserved2;
  unsigned int bfOffBits;
} BITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER {
  unsigned int biSize;
  int biWidth;
  int biHeight;
  unsigned short int biPlanes;
  unsigned short int biBitCount;
  unsigned int biCompression;
  unsigned int biSizeImage;
  int biXPelsPerMeter;
  int biYPelsPerMeter;
  unsigned int biClrUsed;
  unsigned int biClrImportant;
} BITMAPINFOHEADER;

typedef struct tagRGBTRIPLE {
  unsigned char rgbtBlue;
  unsigned char rgbtGreen;
  unsigned char rgbtRed;
} RGBTRIPLE;

typedef struct tagRGBQUAD {
  unsigned char rgbBlue;
  unsigned char rgbGreen;
  unsigned char rgbRed;
  unsigned char rgbReserved;
} RGBQUAD;

#endif

#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L

#define BMP_BF_TYPE 0x4D42 /* word BM */
#define BMP_BF_OFF_BITS 54 /* 14 for file header + 40 for info header (not sizeof(), but packed size) */
#define BMP_BI_SIZE 40 /* packed size of info header */


static unsigned short int WordReadLE(FILE *fp)
{
  // Read a unsigned short int from a file in little endian format 
  unsigned short int lsb, msb;
  lsb = getc(fp);
  msb = getc(fp);
  return (msb << 8) | lsb;
}



static void WordWriteLE(unsigned short int x, FILE *fp)
{
  // Write a unsigned short int to a file in little endian format
  unsigned char lsb = (unsigned char) (x & 0x00FF); putc(lsb, fp); 
  unsigned char msb = (unsigned char) (x >> 8); putc(msb, fp);
}



static unsigned int DWordReadLE(FILE *fp)
{
  // Read a unsigned int word from a file in little endian format 
  unsigned int b1 = getc(fp);
  unsigned int b2 = getc(fp);
  unsigned int b3 = getc(fp);
  unsigned int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void DWordWriteLE(unsigned int x, FILE *fp)
{
  // Write a unsigned int to a file in little endian format 
  unsigned char b1 = (x & 0x000000FF); putc(b1, fp);
  unsigned char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  unsigned char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  unsigned char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



static int LongReadLE(FILE *fp)
{
  // Read a int word from a file in little endian format 
  int b1 = getc(fp);
  int b2 = getc(fp);
  int b3 = getc(fp);
  int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void LongWriteLE(int x, FILE *fp)
{
  // Write a int to a file in little endian format 
  char b1 = (x & 0x000000FF); putc(b1, fp);
  char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



int R2Image::
ReadBMP(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  /* Read file header */
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = WordReadLE(fp);
  bmfh.bfSize = DWordReadLE(fp);
  bmfh.bfReserved1 = WordReadLE(fp);
  bmfh.bfReserved2 = WordReadLE(fp);
  bmfh.bfOffBits = DWordReadLE(fp);
  
  /* Check file header */
  assert(bmfh.bfType == BMP_BF_TYPE);
  /* ignore bmfh.bfSize */
  /* ignore bmfh.bfReserved1 */
  /* ignore bmfh.bfReserved2 */
  assert(bmfh.bfOffBits == BMP_BF_OFF_BITS);
  
  /* Read info header */
  BITMAPINFOHEADER bmih;
  bmih.biSize = DWordReadLE(fp);
  bmih.biWidth = LongReadLE(fp);
  bmih.biHeight = LongReadLE(fp);
  bmih.biPlanes = WordReadLE(fp);
  bmih.biBitCount = WordReadLE(fp);
  bmih.biCompression = DWordReadLE(fp);
  bmih.biSizeImage = DWordReadLE(fp);
  bmih.biXPelsPerMeter = LongReadLE(fp);
  bmih.biYPelsPerMeter = LongReadLE(fp);
  bmih.biClrUsed = DWordReadLE(fp);
  bmih.biClrImportant = DWordReadLE(fp);
  
  // Check info header 
  assert(bmih.biSize == BMP_BI_SIZE);
  assert(bmih.biWidth > 0);
  assert(bmih.biHeight > 0);
  assert(bmih.biPlanes == 1);
  assert(bmih.biBitCount == 24);  /* RGB */
  assert(bmih.biCompression == BI_RGB);   /* RGB */
  int lineLength = bmih.biWidth * 3;  /* RGB */
  if ((lineLength % 4) != 0) lineLength = (lineLength / 4 + 1) * 4;
  assert(bmih.biSizeImage == (unsigned int) lineLength * (unsigned int) bmih.biHeight);

  // Assign width, height, and number of pixels
  width = bmih.biWidth;
  height = bmih.biHeight;
  npixels = width * height;

  // Allocate unsigned char buffer for reading pixels
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = bmih.biSizeImage;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Read buffer 
  fseek(fp, (long) bmfh.bfOffBits, SEEK_SET);
  if (fread(buffer, 1, bmih.biSizeImage, fp) != bmih.biSizeImage) {
    fprintf(stderr, "Error while reading BMP file %s", filename);
    return 0;
  }

  // Close file
  fclose(fp);

  // Allocate pixels for image
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double b = (double) *(p++) / 255;
      double g = (double) *(p++) / 255;
      double r = (double) *(p++) / 255;
      R2Pixel pixel(r, g, b, 1);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
}



int R2Image::
WriteBMP(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Compute number of bytes in row
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;

  // Write file header 
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = BMP_BF_TYPE;
  bmfh.bfSize = BMP_BF_OFF_BITS + rowsize * height;
  bmfh.bfReserved1 = 0;
  bmfh.bfReserved2 = 0;
  bmfh.bfOffBits = BMP_BF_OFF_BITS;
  WordWriteLE(bmfh.bfType, fp);
  DWordWriteLE(bmfh.bfSize, fp);
  WordWriteLE(bmfh.bfReserved1, fp);
  WordWriteLE(bmfh.bfReserved2, fp);
  DWordWriteLE(bmfh.bfOffBits, fp);

  // Write info header 
  BITMAPINFOHEADER bmih;
  bmih.biSize = BMP_BI_SIZE;
  bmih.biWidth = width;
  bmih.biHeight = height;
  bmih.biPlanes = 1;
  bmih.biBitCount = 24;       /* RGB */
  bmih.biCompression = BI_RGB;    /* RGB */
  bmih.biSizeImage = rowsize * (unsigned int) bmih.biHeight;  /* RGB */
  bmih.biXPelsPerMeter = 2925;
  bmih.biYPelsPerMeter = 2925;
  bmih.biClrUsed = 0;
  bmih.biClrImportant = 0;
  DWordWriteLE(bmih.biSize, fp);
  LongWriteLE(bmih.biWidth, fp);
  LongWriteLE(bmih.biHeight, fp);
  WordWriteLE(bmih.biPlanes, fp);
  WordWriteLE(bmih.biBitCount, fp);
  DWordWriteLE(bmih.biCompression, fp);
  DWordWriteLE(bmih.biSizeImage, fp);
  LongWriteLE(bmih.biXPelsPerMeter, fp);
  LongWriteLE(bmih.biYPelsPerMeter, fp);
  DWordWriteLE(bmih.biClrUsed, fp);
  DWordWriteLE(bmih.biClrImportant, fp);

  // Write image, swapping blue and red in each pixel
  int pad = rowsize - width * 3;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      double r = 255.0 * pixel.Red();
      double g = 255.0 * pixel.Green();
      double b = 255.0 * pixel.Blue();
      if (r >= 255) r = 255;
      if (g >= 255) g = 255;
      if (b >= 255) b = 255;
      fputc((unsigned char) b, fp);
      fputc((unsigned char) g, fp);
      fputc((unsigned char) r, fp);
    }

    // Pad row
    for (int i = 0; i < pad; i++) fputc(0, fp);
  }
  
  // Close file
  fclose(fp);

  // Return success
  return 1;  
}



////////////////////////////////////////////////////////////////////////
// PPM I/O
////////////////////////////////////////////////////////////////////////

int R2Image::
ReadPPM(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Read PPM file magic identifier
  char buffer[128];
  if (!fgets(buffer, 128, fp)) {
    fprintf(stderr, "Unable to read magic id in PPM file");
    fclose(fp);
    return 0;
  }

  // skip comments
  int c = getc(fp);
  while (c == '#') {
    while (c != '\n') c = getc(fp);
    c = getc(fp);
  }
  ungetc(c, fp);

  // Read width and height
  if (fscanf(fp, "%d%d", &width, &height) != 2) {
    fprintf(stderr, "Unable to read width and height in PPM file");
    fclose(fp);
    return 0;
  }
  
  // Read max value
  double max_value;
  if (fscanf(fp, "%lf", &max_value) != 1) {
    fprintf(stderr, "Unable to read max_value in PPM file");
    fclose(fp);
    return 0;
  }
  
  // Allocate image pixels
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for PPM file");
    fclose(fp);
    return 0;
  }

  // Check if raw or ascii file
  if (!strcmp(buffer, "P6\n")) {
    // Read up to one character of whitespace (\n) after max_value
    int c = getc(fp);
    if (!isspace(c)) putc(c, fp);

    // Read raw image data 
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
        double r = (double) getc(fp) / max_value;
        double g = (double) getc(fp) / max_value;
        double b = (double) getc(fp) / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }
  else {
    // Read asci image data 
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
  // Read pixel values
  int red, green, blue;
  if (fscanf(fp, "%d%d%d", &red, &green, &blue) != 3) {
    fprintf(stderr, "Unable to read data at (%d,%d) in PPM file", i, j);
    fclose(fp);
    return 0;
  }

  // Assign pixel values
  double r = (double) red / max_value;
  double g = (double) green / max_value;
  double b = (double) blue / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R2Image::
WritePPM(const char *filename, int ascii) const
{
  // Check type
  if (ascii) {
    // Open file
    FILE *fp = fopen(filename, "w");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s", filename);
      return 0;
    }

    // Print PPM image file 
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P3\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%-3d %-3d %-3d  ", r, g, b);
        if (((i+1) % 4) == 0) fprintf(fp, "\n");
      }
      if ((width % 4) != 0) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // Close file
    fclose(fp);
  }
  else {
    // Open file
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s", filename);
      return 0;
    }
    
    // Print PPM image file 
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%c%c%c", r, g, b);
      }
    }
    
    // Close file
    fclose(fp);
  }

  // Return success
  return 1;  
}



////////////////////////////////////////////////////////////////////////
// JPEG I/O
////////////////////////////////////////////////////////////////////////


// #define USE_JPEG
#ifdef USE_JPEG
  extern "C" { 
#   define XMD_H // Otherwise, a conflict with INT32
#   undef FAR // Otherwise, a conflict with windows.h
#   include "jpeg/jpeglib.h"
  };
#endif



int R2Image::
ReadJPEG(const char *filename)
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Initialize decompression info
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo, fp);
  jpeg_read_header(&cinfo, TRUE);
  jpeg_start_decompress(&cinfo);

  // Remember image attributes
  width = cinfo.output_width;
  height = cinfo.output_height;
  npixels = width * height;
  int ncomponents = cinfo.output_components;

  // Allocate pixels for image
  pixels = new R2Pixel [ npixels ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Allocate unsigned char buffer for reading image
  int rowsize = ncomponents * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Read scan lines 
  // First jpeg pixel is top-left, so read pixels in opposite scan-line order
  while (cinfo.output_scanline < cinfo.output_height) {
    int scanline = cinfo.output_height - cinfo.output_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_read_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);

  // Close file
  fclose(fp);

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double r, g, b, a;
      if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
        p++;
      }
      else if (ncomponents == 3) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 4) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = (double) *(p++) / 255;
      }
      else {
        fprintf(stderr, "Unrecognized number of components in jpeg image: %d\n", ncomponents);
        return 0;
      }
      R2Pixel pixel(r, g, b, a);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}


  

int R2Image::
WriteJPEG(const char *filename) const
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Initialize compression info
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, fp);
  cinfo.image_width = width;  /* image width and height, in pixels */
  cinfo.image_height = height;
  cinfo.input_components = 3;   /* # of color components per pixel */
  cinfo.in_color_space = JCS_RGB;   /* colorspace of input image */
  cinfo.dct_method = JDCT_ISLOW;
  jpeg_set_defaults(&cinfo);
  cinfo.optimize_coding = TRUE;
  jpeg_set_quality(&cinfo, 95, TRUE);
  jpeg_start_compress(&cinfo, TRUE);
  
  // Allocate unsigned char buffer for reading image
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Fill buffer with pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      int r = (int) (255 * pixel.Red());
      int g = (int) (255 * pixel.Green());
      int b = (int) (255 * pixel.Blue());
      if (r > 255) r = 255;
      if (g > 255) g = 255;
      if (b > 255) b = 255;
      *(p++) = r;
      *(p++) = g;
      *(p++) = b;
    }
  }



  // Output scan lines
  // First jpeg pixel is top-left, so write in opposite scan-line order
  while (cinfo.next_scanline < cinfo.image_height) {
    int scanline = cinfo.image_height - cinfo.next_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_write_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);

  // Close file
  fclose(fp);

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return number of bytes written
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}




