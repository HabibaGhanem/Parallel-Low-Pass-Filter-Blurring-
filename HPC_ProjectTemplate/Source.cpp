#include <iostream>
#include <math.h>
#include <stdlib.h>
#include<string.h>
#include<msclr\marshal_cppstd.h>
#include <ctime>// include this header 
#include <mpi.h>
#include <cmath>
#pragma once

#using <mscorlib.dll>
#using <System.dll>
#using <System.Drawing.dll>
#using <System.Windows.Forms.dll>
using namespace std;
using namespace msclr::interop;

int* inputImage(int* w, int* h, System::String^ imagePath) //put the size of image in w & h
{
	int* input;


	int OriginalImageWidth, OriginalImageHeight;

	//*********************************************************Read Image and save it to local arrayss*************************	
	//Read Image and save it to local arrayss

	System::Drawing::Bitmap BM(imagePath);

	OriginalImageWidth = BM.Width;
	OriginalImageHeight = BM.Height;
	*w = BM.Width;
	*h = BM.Height;
	int* Red = new int[BM.Height * BM.Width];
	int* Green = new int[BM.Height * BM.Width];
	int* Blue = new int[BM.Height * BM.Width];
	input = new int[BM.Height * BM.Width];
	for (int i = 0; i < BM.Height; i++)
	{
		for (int j = 0; j < BM.Width; j++)
		{
			System::Drawing::Color c = BM.GetPixel(j, i);

			Red[i * BM.Width + j] = c.R;
			Blue[i * BM.Width + j] = c.B;
			Green[i * BM.Width + j] = c.G;

			input[i * BM.Width + j] = ((c.R + c.B + c.G) / 3); //gray scale value equals the average of RGB values

		}

	}
	return input;
}


void createImage(int* image, int width, int height, int index)
{
	System::Drawing::Bitmap MyNewImage(width, height);


	for (int i = 0; i < MyNewImage.Height; i++)
	{
		for (int j = 0; j < MyNewImage.Width; j++)
		{
			//i * OriginalImageWidth + j
			if (image[i * width + j] < 0)
			{
				image[i * width + j] = 0;
			}
			if (image[i * width + j] > 255)
			{
				image[i * width + j] = 255;
			}
			System::Drawing::Color c = System::Drawing::Color::FromArgb(image[i * MyNewImage.Width + j], image[i * MyNewImage.Width + j], image[i * MyNewImage.Width + j]);
			MyNewImage.SetPixel(j, i, c);
		}
	}
	MyNewImage.Save("..//Data//Output//outputRes" + index + ".png");
	//MyNewImage.Save("..//Data//Output//outputRes" + index + ".jpg");
	cout << "result Image Saved " << index << endl;
}


int main()
{
	MPI_Init(NULL, NULL);

	int  size;
	int rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	if (size == 1) {
		int ImageWidth = 4, ImageHeight = 4;

		int start_s, stop_s, TotalTime = 0;

		
		start_s = clock();
		for (int k = 0; k < 10; k++) {
			System::String^ imagePath;
			std::string img;
			img = "..//Data//Input//5N.png";

			imagePath = marshal_as<System::String^>(img);
			int* imageData = inputImage(&ImageWidth, &ImageHeight, imagePath);
			//start_s = clock();
			//SEQUENTIAL CODE

			int rows = ImageHeight + 2;
			int cols = ImageWidth + 2;
			int** image2d = new int* [rows];
			int counter = 0;

			for (int i = 0; i < rows; i++) {
				image2d[i] = new int[cols];
			}
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {
					if (i == 0 || i == (rows - 1) || j == 0 || j == (cols - 1)) {
						image2d[i][j] = 0;
					}
					else {
						image2d[i][j] = imageData[counter];
						counter++;
					}
				}
			}

			int** kernel = new int* [3];
			for (int i = 0; i < 3; i++) {
				kernel[i] = new int[3];
			}
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					kernel[i][j] = 1;
				}
			}

			float** NeWimage2d = new float* [rows];
			for (int i = 0; i < rows; i++)
			{
				NeWimage2d[i] = new float[cols];
			}

			for (int i = 0; i < rows; i++) {

				for (int j = 0; j < cols; j++) {

					NeWimage2d[i][j] = image2d[i][j];
				}
			}

			for (int i = 1; i < rows - 1; i++) {
				for (int j = 1; j < cols - 1; j++) {
					NeWimage2d[i][j] = (image2d[i - 1][j - 1] + image2d[i - 1][j] + image2d[i - 1][j + 1]
						+ image2d[i][j - 1] + image2d[i][j] + image2d[i][j + 1]
						+ image2d[i + 1][j - 1] + image2d[i + 1][j] + image2d[i + 1][j + 1]) / 9.0;
					NeWimage2d[i][j] = round(NeWimage2d[i][j]);

				}

			}

			for (int i = 1; i < rows - 1; i++) {
				for (int j = 1; j < cols - 1; j++) {
					imageData[(i - 1) * (cols - 2) + (j - 1)] = NeWimage2d[i][j];
				}
			}

			createImage(imageData, ImageWidth, ImageHeight, k);
			//stop_s = clock();
			free(imageData);
		}
		stop_s = clock();
		TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;
		cout << "time: " << TotalTime << endl;

		//free(imageData);

	}
	else {
		int ImageWidth, ImageHeight;

		int start_s, stop_s, TotalTime = 0;

		int rows, cols; //height and width plus padding

		int reminder;//reminder of rows after equally distributed them on size-1 ranks "will be computed by last rank"
		int SizePerRank; //# of elements that will be computed by each rank (last rank is excluded)
		int RecvCount;//total # of elements send to each rank ,, no of elements needed by each processor 	
		int* LocalData = NULL;
		float* LocalDataSum = NULL;//LocalData of excatly # of rows and cols (ImageHeight & ImageWidth)after computed
		int TotalRowsPerRank; //Total # of rows send to each rank

		int RowsPerRank;//# of rows to computed send to each rank
		int* RecvImage = NULL;

		MPI_Status stat;
		start_s = clock();
		for (int k = 0; k < 10; k++) {
			if (rank == 0) {
				ImageWidth = 4;
				ImageHeight = 4;


				System::String^ imagePath;
				std::string img;
				img = "..//Data//Input//5N.png";

				imagePath = marshal_as<System::String^>(img);
				int* imageData = inputImage(&ImageWidth, &ImageHeight, imagePath);

				//start_s = clock();

				rows = ImageHeight + 2;
				cols = ImageWidth + 2;

				int** image2d = new int* [rows];

				for (int i = 0; i < rows; i++) {
					image2d[i] = new int[cols];
				}
				int counter = 0;
				for (int i = 0; i < rows; i++)
				{
					for (int j = 0; j < cols; j++)
					{
						if (i == 0 || i == (rows - 1) || j == 0 || j == (cols - 1))
						{
							image2d[i][j] = 0;
						}
						else
						{
							image2d[i][j] = imageData[counter];
							counter++;
						}
					}
				}
				/*
				cout << "image1d size : " << ImageHeight * ImageWidth << endl;
				for (int i = 0; i < ImageHeight * ImageWidth; i++) {
					imageData[i];
					cout << "image1d pixel : " << i<< "  " << imageData[i] << endl;
				}
				*/
				//cout << "After padding" << endl;
				int* image1d = new int[rows * cols]; //padding + image1d are right
				//cout << "image1d size : " << rows * cols << endl;
				for (int i = 0; i < rows; i++) {
					for (int j = 0; j < cols; j++) {
						image1d[i * cols + j] = image2d[i][j];
						//cout << "image1d pixel : " << i * cols + j << "  " << image1d[i * cols + j] << endl;
					}
				}
				reminder = (ImageHeight % (size - 1)) * cols;
				SizePerRank = (ImageHeight / (size - 1)) * cols;

				for (int i = 1; i < size; i++)
				{
					if (i == size - 1)
					{
						RecvCount = reminder + cols * 2;
						MPI_Send(&reminder, 1, MPI_INT, i, i, MPI_COMM_WORLD);
					}
					else
					{
						RecvCount = SizePerRank + cols * 2;
						MPI_Send(&SizePerRank, 1, MPI_INT, i, i, MPI_COMM_WORLD);
					}
					MPI_Send(&RecvCount, 1, MPI_INT, i, i, MPI_COMM_WORLD); //Send and received successfully
					//cout << rank << " Rank   ";
					//cout << "RecvCount : " << RecvCount << endl;							
					MPI_Send(&image1d[(SizePerRank * i)], RecvCount, MPI_INT, i, i, MPI_COMM_WORLD);
					//cout << "start index of rank : " << i << "  is  " << SizePerRank * i; //start index of each rank is right
					MPI_Send(&cols, 1, MPI_INT, i, i, MPI_COMM_WORLD);
					MPI_Send(&ImageWidth, 1, MPI_INT, i, i, MPI_COMM_WORLD);
				}

				RecvCount = SizePerRank + cols * 2;
				//cout << "Rank : " << rank << "   RecvCount : " << RecvCount; //right

				TotalRowsPerRank = RecvCount / cols;
				RowsPerRank = SizePerRank / cols;
				LocalDataSum = new float[RowsPerRank * ImageWidth];

				int Counter = 0;
				for (int i = 1; i < TotalRowsPerRank - 1; i++) { //calculated successfully and in the right positions
					for (int j = 1; j < cols - 1; j++)
					{
						LocalDataSum[Counter] = (image1d[(i - 1) * cols + (j - 1)] + image1d[(i - 1) * cols + j]
							+ image1d[(i - 1) * cols + (j + 1)] + image1d[i * cols + (j - 1)]
							+ image1d[i * cols + j] + image1d[i * cols + (j + 1)]
							+ image1d[(i + 1) * cols + (j - 1)] + image1d[(i + 1) * cols + j]
							+ image1d[(i + 1) * cols + (j + 1)]) / 9.0;
						LocalDataSum[Counter] = round(LocalDataSum[Counter]);
						Counter++;
					}
				}
				/*
				cout << rank << "  Rank" << endl;
				for (int i = 0; i < RowsPerRank * ImageWidth; i++) {
					cout << LocalDataSum[i] << " , ";
				}
				cout << endl;
				*/
				int inc = 0;
				RecvImage = new int[ImageHeight * ImageWidth];
				for (int i = 0; i < RowsPerRank * ImageWidth; i++) {
					RecvImage[inc++] = LocalDataSum[i];
				}

				int s;
				for (int j = 1; j < size; j++) {
					MPI_Recv(&s, 1, MPI_INT, j, j, MPI_COMM_WORLD, &stat);
					if (s == 0)continue;
					MPI_Recv(LocalDataSum, s, MPI_INT, j, j, MPI_COMM_WORLD, &stat);
					for (int i = 0; i < s; i++)
						RecvImage[inc++] = LocalDataSum[i];
				}
				/*
				for (int i = 0; i < ImageWidth * ImageHeight; i++) {
					cout << RecvImage[i] << " , ";
				}
				*/

				createImage(RecvImage, ImageWidth, ImageHeight, k);
				//stop_s = clock();
				//TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;
				//cout << "time: " << TotalTime << endl;
				free(imageData);

			}
			else {
				MPI_Recv(&SizePerRank, 1, MPI_INT, 0, rank, MPI_COMM_WORLD, &stat);
				if (SizePerRank != 0) {
					MPI_Recv(&RecvCount, 1, MPI_INT, 0, rank, MPI_COMM_WORLD, &stat);
					//cout << "Rank : " << rank << "   RecvCount : " << RecvCount << endl;

					LocalData = new int[RecvCount];
					MPI_Recv(LocalData, RecvCount, MPI_INT, 0, rank, MPI_COMM_WORLD, &stat);
					/* //received data successfully
					for (int i = 0; i < RecvCount; i++) {
						if(i==0)cout << "LocalData of rank : " << rank  <<" are ";
						cout << LocalData[i] << " , ";
					}*/
					MPI_Recv(&cols, 1, MPI_INT, 0, rank, MPI_COMM_WORLD, &stat);
					MPI_Recv(&ImageWidth, 1, MPI_INT, 0, rank, MPI_COMM_WORLD, &stat);


					TotalRowsPerRank = RecvCount / cols;
					RowsPerRank = SizePerRank / cols;
					LocalDataSum = new float[RowsPerRank * ImageWidth];

					int Counter = 0;
					for (int i = 1; i < TotalRowsPerRank - 1; i++) { //calculated successfully and in the right positions
						for (int j = 1; j < cols - 1; j++)
						{
							LocalDataSum[Counter] = (LocalData[(i - 1) * cols + (j - 1)] + LocalData[(i - 1) * cols + j]
								+ LocalData[(i - 1) * cols + (j + 1)] + LocalData[i * cols + (j - 1)]
								+ LocalData[i * cols + j] + LocalData[i * cols + (j + 1)]
								+ LocalData[(i + 1) * cols + (j - 1)] + LocalData[(i + 1) * cols + j]
								+ LocalData[(i + 1) * cols + (j + 1)]) / 9.0;
							LocalDataSum[Counter] = round(LocalDataSum[Counter]);
							Counter++;
						}
					}
					/*
					cout << rank << "  Rank" << endl;
					for (int i = 0; i < RowsPerRank * ImageWidth; i++) {
						cout << LocalDataSum[i] << " , ";
					}
					cout << endl;
					*/
					int s = RowsPerRank * ImageWidth; //exactly # of cells computed by each rank
					MPI_Send(&s, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
					MPI_Send(LocalDataSum, RowsPerRank * ImageWidth, MPI_INT, 0, rank, MPI_COMM_WORLD);
				}
				else {
					int s = 0;
					MPI_Send(&s, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
				}
			}
		}
		if (rank == 0) {
			stop_s = clock();
			TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;
			cout << "time: " << TotalTime << endl;
		}
	}
	MPI_Finalize();
	return 0;
}



