#include <opencv2/opencv.hpp>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <algorithm>
#include <math.h>
using namespace cv;
using namespace std;

int EstimateIlluminantGrey(Mat I, double p) {

	int Ic = 0, L = 256;
	int width = I.cols, height = I.rows;
	
	int pixelTh = 0;
	int histI[256] = { 0 };

	pixelTh = int(p * width*height);
	//histI=imhist(I);
	//calc image histogram,I =B or G or R channel
	for (int i = 0; i < height; i++){
		for (int j = 0; j < width; j++){
			int pixel = I.at<uchar>(i, j);
			histI[pixel]++;
		}
	}
	//get max&min value 	
	//minMaxIdx(I, &minVal, &maxVal);
	int minVal = I.at<uchar>(0, 0), maxVal = I.at<uchar>(0, 0);
	 for (int i = 0; i < height; i++){
		for (int j = 0; j < width; j++){
			if(I.at<uchar>(i, j)> maxVal){
				maxVal= I.at<uchar>(i, j);
			}
			if(I.at<uchar>(i, j)< minVal){
				minVal= I.at<uchar>(i, j);
			}
			
		}
	 }	
	int sum = 0;
	for (int k = maxVal; k >= minVal; k--)
	{

		sum = sum + histI[k];

		if (sum > pixelTh)
		{
			Ic = k;
			break;
		}

	}
	return Ic;
}

double EstimateCCT(int R, int G, int B) {

	//Constant parameters from 
	double A0 = -949.86315;
	double A1 = 6253.80338;
	double A2 = 28.70599;
	double A3 = 0.00004;

	double t1 = 0.92159;
	double t2 = 0.20039;
	double t3 = 0.07125;

	double xe = 0.3366;
	double ye = 0.1735;
	double x = 0, y = 0;
	double H = 0, CCT = 0;
	//Calculate x and y from estimated illuminant values
	/*
	XYZ_Conv_matrix = [ 0.4124 0.3576 0.1805;
						0.2126 0.7152 0.0722;
						0.0193 0.152 0.9505];
	XYZ = XYZ_Conv_matrix * double(iEstm');
	*/
	double X = 0.4124*R + 0.3576*G + 0.1805*B;
	double Y = 0.2126*R + 0.7152*G + 0.0722*B;
	double Z = 0.0193*R + 0.152*G + 0.9505*B;
	x = X / (X + Y + Z);
	y = Y / (X + Y + Z);

	H = -((x - xe) / (y - ye));

	CCT = A0 + (A1*exp(H / t1)) + (A2*exp(H / t2)) + (A3*exp(H / t3));
	
	return CCT;

}
Mat performAWB(Mat src, double gray_percent) {

	Mat BGRChannel[3];
	//Mat dst;
	int width = src.cols, height = src.rows;
	Mat dst(height, width, CV_8UC3);
	split(src, BGRChannel);
	int Bc, Gc, Rc;
	Bc = EstimateIlluminantGrey(BGRChannel[0], gray_percent);
	Gc = EstimateIlluminantGrey(BGRChannel[1], gray_percent);
	Rc = EstimateIlluminantGrey(BGRChannel[2], gray_percent);
	cout << "Bc  " << Bc << endl;
	cout << "Gc  " << Gc << endl;
	cout << "Rc  " << Rc << endl;
	//calc current_CCT and reference_CCT
	double CCT_Estm = 0, CCT_Ref = 0;
	CCT_Estm = EstimateCCT(Rc, Gc, Bc);
	cout << "CT_estimate " << CCT_Estm <<"K"<< endl;
	CCT_Ref = EstimateCCT(Gc, Gc, Gc);
	cout << "CT_reference " << CCT_Ref <<"K"<< endl;
	//int iEstm[] = {Rc,Gc,Bc};
	//int iRef[] = {Gc,Gc,Gc};
	double Kr, Kg, Kb;
	Kr = (double)Gc / Rc;//ref_R/estm_R
	Kg = (double)Gc / Gc;//ref_G/estm_G
	Kb = (double)Gc / Bc;//ref_B/estm_B

	cout << "Kr " << Kr << endl;
	cout << "Kg " << Kg << endl;
	cout << "Kb " << Kb << endl;

	//Mat K = (Mat_<double>(3,3)<< Kr, 0, 0, 0, Kg, 0, 0, 0, Kb);
	//Mat K = (Mat_<double>(3, 3) << Kb, 0, 0, 0, Kg, 0, 0, 0, Kr);

	double Tr, Tg, Tb;
	Tr = fmax(1, (CCT_Estm - CCT_Ref) / 100) * (Kr - 1);
	Tg = 0;
	Tb = fmax(1, (CCT_Ref - CCT_Estm) / 100) * (Kb - 1);
	cout << "Tr " << Tr << endl;
	cout << "Tg " << Tg << endl;
	cout << "Tb " << Tb << endl;
	//Mat T = (Mat_<double>(3,1)<< Tr, Tg, Tb);
	//Mat T = (Mat_<double>(3, 1) << Tb, Tg, Tr);
    int FWB[3];
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			//Vec<int,3>  Fxy(src.at<Vec3b>(i, j)[0], src.at<Vec3b>(i, j)[1], src.at<Vec3b>(i, j)[2])£»
			//Mat Fxy = (Mat_<double>(3, 1) << src.at<Vec3b>(i, j)[0], src.at<Vec3b>(i, j)[1], src.at<Vec3b>(i, j)[2]);
			//Mat FWB = int(K * Fxy + T);
			FWB[0] = round(Kb * src.at<Vec3b>(i, j)[0]+Tb);
			FWB[1] = round(Kg * src.at<Vec3b>(i, j)[1]+Tg);
			FWB[2] = round(Kr * src.at<Vec3b>(i, j)[2]+Tr);
			
			for (int n = 0; n < 3; n++) {
				
				if (FWB[n] > 255) {
					FWB[n] = 255;
				}
				else if (FWB[n] < 0) {
					FWB[n] = 0;
				}
			
				dst.at<Vec3b>(i, j)[n] = FWB[n];
				/*
				if (dst.at<Vec3b>(i,j)[n] > 255) {
					dst.at<Vec3b>(i,j)[n] = 255;
				}else if (dst.at<Vec3b>(i, j)[n] < 0){
					dst.at<Vec3b>(i, j)[n] = 0;
				}
				*/
			}
		}
	}

	return dst;

}


int main() {

	Mat src = imread("C:\\Users\\ggcle\\Desktop\\AWB\\AWB Code\\image4.bmp");
	if (src.empty())
	{
		cout << "Could not open src " << endl;
		return -1;
	}
	Mat dst = performAWB(src, 0.1);
	//imshow("origin", src);
	//imshow("result", dst);
	//for display two images in single window//////
	Mat Comp(Size(src.cols * 2, src.rows), src.type(), Scalar::all(0));
	Mat Roi = Comp(Rect(0, 0, src.cols, src.rows));
	src.copyTo(Roi);
	Roi = Comp(Rect(src.cols, 0, src.cols, src.rows));
	dst.copyTo(Roi);
	//imshow("original vs result", Comp);
	imwrite("result4.jpg", Comp);
	/////////////
	//imwrite("../result.jpg", dst);
	waitKey(0);
	destroyAllWindows();


	return 0;
}


