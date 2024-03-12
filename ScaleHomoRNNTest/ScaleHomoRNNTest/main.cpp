#include "pch.h"
#include "seal/seal.h"
#include <fstream>
#include <iomanip>
#include <chrono>

using namespace std;
using namespace chrono;
using namespace seal;

#define row 8000
#define col 39
#define time 48
#define hidenode 31
#define slotCount 4096
#define partyNum 3
#define batchsize 250
#define levelZeroScale pow(2, 31)
#define levelOne { 16588149047881338678, 5073713403208523954, 11766700155144244509, 6013047097330863554 }
#define levelOneScale 2147762211.1294265
#define levelTwo { 13212098252820938193, 1095622039496364029, 6909526250636377009, 14736031212968045950 }
#define levelTwoScale 2148450593.8227954
#define levelThree { 18295512784679830468, 15921918161301782861, 17525725753082260777, 13618162003755478695 }
#define levelThreeScale 2149942862.4245496
#define levelFour { 1764438775781066231, 13275548663487857381, 12258101589568993885, 6939347545201792930 }
#define levelFourScale 2153176981.5671215

int mod(int x, int y)
{
	return x % y < 0 ? x % y + y : x % y;
}

void parameterEncode(Ciphertext& parameter, int rowN, int colN, int squN, Evaluator& evaluator, CKKSEncoder& ckks_encoder, vector<double> empty, double w[][hidenode], double b[] = new double{ 1000000 })
{
	int finRow = ceil(rowN * 1.0 / squN) * squN;
	int finCol = ceil(colN * 1.0 / squN) * squN;
	int copyNum = log2(slotCount / ((long long)finRow * finCol));
	vector<double> temp = empty;
	Plaintext tempPlain;
	for (int i = 0; i < (finRow * finCol) / (squN * squN); i++)
	{
		for (int j = 0; j < squN; j++)
		{
			for (int k = 0; k < squN; k++)
			{
				if ((i / (finCol / squN)) * squN + j < rowN)
				{
					if ((i % (finCol / squN)) * squN + k < colN)
					{
						for (int m = 0; m < pow(2, copyNum); m++)
						{
							temp[(long long)j * (slotCount / squN) + mod(k - j, squN) + (i * pow(2, copyNum) + m) * squN] = w[(i % (finCol / squN)) * squN + k][(i / (finCol / squN)) * squN + j] * 3.735111370627234931;
						}
					}
					else if ((i % (finCol / squN)) * squN + k == colN)
					{
						if (b[0] != 1000000)
						{
							for (int m = 0; m < pow(2, copyNum); m++)
							{
								temp[(long long)j * (slotCount / squN) + mod(k - j, squN) + (i * pow(2, copyNum) + m) * squN] = b[(i / (finCol / squN)) * squN + j] * 3.735111370627234931;
							}
						}
					}
				}
			}
		}
	}
	ckks_encoder.encode(temp, levelZeroScale, tempPlain);
	evaluator.add_plain_inplace(parameter, tempPlain);
}

void parameterEncodeTranspose(Ciphertext& parameter, int rowN, int colN, int squN, Evaluator& evaluator, CKKSEncoder& ckks_encoder, vector<double> empty, double w[][hidenode], double b[] = new double{ 1000000 })
{
	int finRow = ceil(rowN * 1.0 / squN) * squN;
	int finCol = ceil(colN * 1.0 / squN) * squN;
	int copyNum = log2(slotCount / ((long long)finRow * finCol));
	vector<double> temp = empty;
	Plaintext tempPlain;
	for (int i = 0; i < (finRow * finCol) / (squN * squN); i++)
	{
		for (int j = 0; j < squN; j++)
		{
			for (int k = 0; k < squN; k++)
			{
				if ((i / (finCol / squN)) * squN + j < rowN)
				{
					if ((i % (finCol / squN)) * squN + k < colN)
					{
						for (int m = 0; m < pow(2, copyNum); m++)
						{
							temp[(long long)j * (slotCount / squN) + mod(k - j, squN) + (i * pow(2, copyNum) + m) * squN] = w[(i / (finCol / squN)) * squN + j][(i % (finCol / squN)) * squN + k] * 5.282245157314995551;
						}
					}
				}
			}
		}
	}
	ckks_encoder.encode(temp, levelZeroScale, tempPlain);
	evaluator.add_plain_inplace(parameter, tempPlain);
}

void outParameterEncode(Ciphertext& parameter, int rowN, int squN, Evaluator& evaluator, CKKSEncoder& ckks_encoder, vector<double> empty, double w[], double b)
{
	vector<double> temp = empty;
	Plaintext tempPlain;
	for (int i = 0; i < rowN + 1; i++)
	{
		for (int j = 0; j < (slotCount / squN - squN) / 2; j++)
		{
			if (i != rowN)
			{
				temp[(long long)i % squN * slotCount / squN + i / squN * slotCount / (squN * 2) + j] = w[i];
			}
			else
			{
				temp[(long long)i % squN * slotCount / squN + i / squN * slotCount / (squN * 2) + j] = b;
			}
		}
	}
	ckks_encoder.encode(temp, levelZeroScale, tempPlain);
	evaluator.add_plain_inplace(parameter, tempPlain);
}

void parameterRotate(Ciphertext c1[], int rowN, int colN, int squN, Ciphertext result[], Evaluator& evaluator, GaloisKeys gal_key, Plaintext one_Plain)
{
	int finRow = ceil(rowN * 1.0 / squN);
	int finCol = ceil(colN * 1.0 / squN);
	for (int i = 0; i < finRow * finCol; i++)
	{
		result[i * squN] = c1[i];
		if (c1[i].parms_id() == one_Plain.parms_id())
		{
			evaluator.multiply_plain_inplace(result[i * squN], one_Plain);
			evaluator.rescale_to_next_inplace(result[i * squN]);
		}
		for (int j = 1; j < squN; j++)
		{
			evaluator.rotate_vector(result[i * squN + j - 1], 1, gal_key, result[i * squN + j]);
		}
	}
}

void dataEncode(vector<double> data[], int rowN, int colN, int squN, double input[row][col][time], int ifodd, int index)
{
	int finRow = ceil(rowN * 1.0 / squN) * squN;
	int finCol = ceil(colN * 1.0 / squN) * squN;
	int rowIndex = 0;
	int vectorPos = 0;
	for (int i = batchsize * (index % (row / batchsize)); i < batchsize * ((index % (row / batchsize)) + 1); i++)
	{
		for (int j = 0; j < time / 2; j++)
		{
			for (int k = 0; k < finCol; k++)
			{
				if (k < colN)
				{
					data[vectorPos + k / squN + (rowIndex % finRow / squN) * (finCol / squN)][(long long)mod(k - rowIndex, squN) * (slotCount / squN) + (long long)rowIndex / finRow * squN + rowIndex % squN] = input[i][k][2 * j + ifodd] * 3.735111370627234931;
				}
				else if (k == colN)
				{
					data[vectorPos + k / squN + (rowIndex % finRow / squN) * (finCol / squN)][(long long)mod(k - rowIndex, squN) * (slotCount / squN) + (long long)rowIndex / finRow * squN + rowIndex % squN] = 3.735111370627234931;
				}
			}
			rowIndex++;
			if (rowIndex == slotCount * (finRow / squN) * (finCol / squN) / finCol * (slotCount / ((long long)squN * squN) - 1) / (slotCount / ((long long)squN * squN)))
			{
				rowIndex = 0;
				vectorPos += (finRow / squN) * (finCol / squN);
			}
		}
	}
}

void dataCode(vector<double> data[], int rowN, int colN, int squN, double input[row][col][time], int ifodd, int index)
{
	int finRow = ceil(rowN * 1.0 / squN) * squN;
	int finCol = ceil(colN * 1.0 / squN) * squN;
	int rowIndex = 0;
	int vectorPos = 0;
	for (int i = batchsize * (index % (row / batchsize)); i < batchsize * ((index % (row / batchsize)) + 1); i++)
	{
		for (int j = 0; j < time / 2; j++)
		{
			for (int k = 0; k < finCol; k++)
			{
				if (k < colN)
				{
					data[vectorPos + k / squN * (finRow / squN) + (rowIndex % finRow / squN)][k % squN * (data[0].size() / squN) + (long long)rowIndex / finRow * squN + rowIndex % squN] = input[i][k][2 * j + ifodd] * pow(2, 7.75);
				}
				else if (k == colN)
				{
					data[vectorPos + k / squN * (finRow / squN) + (rowIndex % finRow / squN)][k % squN * (data[0].size() / squN) + (long long)rowIndex / finRow * squN + rowIndex % squN] = 1 * pow(2, 7.75);
				}
				else
				{
					data[vectorPos + k / squN * (finRow / squN) + (rowIndex % finRow / squN)][k % squN * (data[0].size() / squN) + (long long)rowIndex / finRow * squN + rowIndex % squN] = 0;
				}
			}
			rowIndex++;
			if (rowIndex == data[0].size() * (finRow / squN) * (finCol / squN) / finCol * (data[0].size() / ((long long)squN * squN) - 1) / (data[0].size() / ((long long)squN * squN)))
			{
				rowIndex = 0;
				vectorPos += (finRow / squN) * (finCol / squN);
			}
		}
	}
}

void matrixMultiple_CipPlain(Ciphertext c1[], vector<double> v2[], int rowN, int colN, int squN, int resultNum, Ciphertext result[], Evaluator& evaluator, CKKSEncoder& ckks_encoder)
{
	int finRow = ceil(rowN * 1.0 / squN);
	int finCol = ceil(colN * 1.0 / squN);
	int inputLen = resultNum / (finRow * finRow);
	Ciphertext tempCip;
	Plaintext** p2 = new Plaintext * [(long long)finRow * finCol * inputLen];
	for (int i = 0; i < finRow * finCol * inputLen; i++)
	{
		p2[i] = new Plaintext[squN];
		ckks_encoder.encode(v2[i], c1[0].parms_id(), c1[0].scale(), p2[i][0]);
		for (int j = 1; j < squN; j++)
		{
			vector<double> temp;
			rotate_copy(v2[i].begin(), v2[i].begin() + (long long)slotCount * j / squN, v2[i].end(), back_inserter(temp));
			ckks_encoder.encode(temp, c1[0].parms_id(), c1[0].scale(), p2[i][j]);
		}
	}
	for (int i = 0; i < inputLen; i++)
	{
		for (int j = 0; j < finRow; j++)
		{
			for (int k = 0; k < finRow; k++)
			{
				for (int m = 0; m < finCol; m++)
				{
					for (int n = 0; n < squN; n++)
					{
						if (m == 0 && n == 0)
						{
							evaluator.multiply_plain(c1[(j * finCol + m) * squN + n], p2[i * finRow * finCol + k * finCol + m][n], result[i * finRow * finRow + j * finRow + k]);
						}
						else
						{
							evaluator.multiply_plain(c1[(j * finCol + m) * squN + n], p2[i * finRow * finCol + k * finCol + m][n], tempCip);
							evaluator.add_inplace(result[i * finRow * finRow + j * finRow + k], tempCip);
						}
					}
				}
			}
		}
	}
	for (int i = 0; i < finRow * finCol * inputLen; i++)
	{
		delete[] p2[i];
	}
	delete[] p2;
}

void matrixMultiple_CipCip(Ciphertext c1[], Ciphertext c2[], int rowN, int colN, int squN, int resultNum, Ciphertext result[], Evaluator& evaluator, GaloisKeys gal_key)
{
	int finRow = ceil(rowN * 1.0 / squN);
	int finCol = ceil(colN * 1.0 / squN);
	int inputLen = resultNum / (finRow * finRow);
	Ciphertext tempCip;
	Ciphertext** cTrans = new Ciphertext * [(long long)finRow * finCol * inputLen];
	for (int i = 0; i < finRow * finCol * inputLen; i++)
	{
		cTrans[i] = new Ciphertext[squN];
		cTrans[i][0] = c2[i];
		for (int j = 1; j < squN; j++)
		{
			evaluator.rotate_vector(cTrans[i][j - 1], slotCount / squN, gal_key, cTrans[i][j]);
		}
	}
	for (int i = 0; i < inputLen; i++)
	{
		for (int j = 0; j < finRow; j++)
		{
			for (int k = 0; k < finRow; k++)
			{
				for (int m = 0; m < finCol; m++)
				{
					for (int n = 0; n < squN; n++)
					{
						if (m == 0 && n == 0)
						{
							evaluator.multiply(c1[(j * finCol + m) * squN + n], cTrans[i * finRow * finCol + k * finCol + m][n], result[i * finRow * finRow + j * finRow + k]);
						}
						else
						{
							evaluator.multiply(c1[(j * finCol + m) * squN + n], cTrans[i * finRow * finCol + k * finCol + m][n], tempCip);
							evaluator.add_inplace(result[i * finRow * finRow + j * finRow + k], tempCip);
						}
					}
				}
			}
		}
	}
	for (int i = 0; i < finRow * finCol * inputLen; i++)
	{
		delete[] cTrans[i];
	}
	delete[] cTrans;
}

void matrixRefresh(vector<double> v1[], vector<double> v2[], int rowN, int colN, int squN, int resultNum, int dataNum)
{
	int finRow = ceil(rowN * 1.0 / squN);
	int finCol = ceil(colN * 1.0 / squN);
	int inputLen = resultNum / (finRow * finCol);
	int colIndex = 0;
	int v1Index;
	int colPos;
	for (int i = 0; i < inputLen; i++)
	{
		for (int j = 0; j < finRow; j++)
		{
			colPos = colN;
			for (int k = 0; k < finCol; k++)
			{
				for (int m = 0; m < v1[0].size(); m++)
				{
					v1Index = (m + colIndex * (v1[0].size() / squN)) % v1[0].size();
					if (v1Index / (v1[0].size() / squN) < colPos)
					{
						v2[i * finRow * finCol + j * finCol + k][m] = v1[i * finRow * finCol + j + k * finCol][v1Index];
					}
					else
					{
						v2[i * finRow * finCol + j * finCol + k][m] = 0;
					}
					colIndex++;
					if (colIndex == squN)
					{
						colIndex = 0;
					}
				}
				colPos -= squN;
			}
		}
	}
}

void dataDivideEncode(vector<double> input[], vector<double> output1[], vector<double> output2[], int dataSize, int rowN, int colN, int inSqu, int outSqu, vector<double> empty, int resultNum, double b_default)
{
	double*** temp = new double** [2];
	for (int i = 0; i < 2; i++)
	{
		temp[i] = new double* [dataSize / 2];
		for (int j = 0; j < dataSize / 2; j++)
		{
			temp[i][j] = new double[colN];
		}
	}
	int inRow = ceil(rowN * 1.0 / inSqu);
	int inCol = ceil(colN * 1.0 / inSqu);
	int rowPos = slotCount * inRow * inCol * (slotCount / (inSqu * inSqu) - 1) / (slotCount / (inSqu * inSqu)) / (inCol * inSqu);
	for (int i = 0; i < dataSize; i++)
	{
		for (int j = 0; j < colN; j++)
		{
			temp[i % 2][i / 2][j] = input[(i / rowPos) * inRow * inCol + (j / inSqu) * inCol + i % (inCol * inSqu) / inSqu][(long long)(i % rowPos / (inCol * inSqu)) * inSqu + i % inSqu + j % inSqu * slotCount / inSqu];
		}
	}
	for (int i = 0; i < resultNum; i++)
	{
		output1[i] = empty;
		output2[i] = empty;
	}
	int outRow = ceil(rowN * 1.0 / outSqu) * outSqu;
	int outCol = ceil(colN * 1.0 / outSqu) * outSqu;
	int rowIndex = 0;
	int vectorPos = 0;
	for (int i = 0; i < dataSize / 2; i++)
	{
		for (int j = 0; j < outCol; j++)
		{
			if (j < colN)
			{
				output1[vectorPos + j / outSqu + (rowIndex % outRow / outSqu) * (outCol / outSqu)][mod(j - rowIndex, outSqu) * (output1[0].size() / outSqu) + (long long)rowIndex / outRow * outSqu + rowIndex % outSqu] = temp[0][i][j] * 3.735111370627234931 / pow(2, 15.5);
				output2[vectorPos + j / outSqu + (rowIndex % outRow / outSqu) * (outCol / outSqu)][mod(j - rowIndex, outSqu) * (output2[0].size() / outSqu) + (long long)rowIndex / outRow * outSqu + rowIndex % outSqu] = temp[1][i][j] * 3.735111370627234931 / pow(2, 15.5);
			}
			else
			{
				output1[vectorPos + j / outSqu + (rowIndex % outRow / outSqu) * (outCol / outSqu)][mod(j - rowIndex, outSqu) * (output1[0].size() / outSqu) + (long long)rowIndex / outRow * outSqu + rowIndex % outSqu] = b_default * 3.735111370627234931;
				output2[vectorPos + j / outSqu + (rowIndex % outRow / outSqu) * (outCol / outSqu)][mod(j - rowIndex, outSqu) * (output2[0].size() / outSqu) + (long long)rowIndex / outRow * outSqu + rowIndex % outSqu] = b_default * 3.735111370627234931;
			}
		}
		rowIndex++;
		if (rowIndex == output1[0].size() * (outRow / outSqu) * (outCol / outSqu) / outCol * (output1[0].size() / ((long long)outSqu * outSqu) - 1) / (output1[0].size() / ((long long)outSqu * outSqu)))
		{
			rowIndex = 0;
			vectorPos += (outRow / outSqu) * (outCol / outSqu);
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < dataSize / 2; j++)
		{
			delete[] temp[i][j];
		}
		delete[] temp[i];
	}
	delete[] temp;
}

void dataDivideEncodeThree(vector<double> input[], vector<double> output1[], vector<double> output2[], vector<double> output3[], int dataSize, int rowN, int colN, int inSqu, int outSqu, vector<double> empty, int resultNum, double b_default)
{
	double*** temp = new double** [3];
	for (int i = 0; i < 3; i++)
	{
		temp[i] = new double* [dataSize / 3];
		for (int j = 0; j < dataSize / 3; j++)
		{
			temp[i][j] = new double[colN];
		}
	}
	int inRow = ceil(rowN * 1.0 / inSqu);
	int inCol = ceil(colN * 1.0 / inSqu);
	int rowPos = slotCount * inRow * inCol * (slotCount / (inSqu * inSqu) - 1) / (slotCount / (inSqu * inSqu)) / (inCol * inSqu);
	for (int i = 0; i < dataSize; i++)
	{
		for (int j = 0; j < colN; j++)
		{
			temp[i % 3][i / 3][j] = input[(i / rowPos) * inRow * inCol + (j / inSqu) * inCol + i % (inCol * inSqu) / inSqu][(long long)(i % rowPos / (inCol * inSqu)) * inSqu + i % inSqu + j % inSqu * slotCount / inSqu];
		}
	}
	for (int i = 0; i < resultNum; i++)
	{
		output1[i] = empty;
		output2[i] = empty;
		output3[i] = empty;
	}
	int outRow = ceil(rowN * 1.0 / outSqu) * outSqu;
	int outCol = ceil(colN * 1.0 / outSqu) * outSqu;
	int rowIndex = 0;
	int vectorPos = 0;
	for (int i = 0; i < dataSize / 3; i++)
	{
		for (int j = 0; j < outCol; j++)
		{
			if (j < colN)
			{
				output1[vectorPos + j / outSqu + (rowIndex % outRow / outSqu) * (outCol / outSqu)][mod(j - rowIndex, outSqu) * (output1[0].size() / outSqu) + (long long)rowIndex / outRow * outSqu + rowIndex % outSqu] = temp[0][i][j] * 3.735111370627234931 / pow(2, 15.5);
				output2[vectorPos + j / outSqu + (rowIndex % outRow / outSqu) * (outCol / outSqu)][mod(j - rowIndex, outSqu) * (output2[0].size() / outSqu) + (long long)rowIndex / outRow * outSqu + rowIndex % outSqu] = temp[1][i][j] * 3.735111370627234931 / pow(2, 15.5);
				output3[vectorPos + j / outSqu + (rowIndex % outRow / outSqu) * (outCol / outSqu)][mod(j - rowIndex, outSqu) * (output3[0].size() / outSqu) + (long long)rowIndex / outRow * outSqu + rowIndex % outSqu] = temp[2][i][j] * 3.735111370627234931 / pow(2, 15.5);
			}
			else
			{
				output1[vectorPos + j / outSqu + (rowIndex % outRow / outSqu) * (outCol / outSqu)][mod(j - rowIndex, outSqu) * (output1[0].size() / outSqu) + (long long)rowIndex / outRow * outSqu + rowIndex % outSqu] = b_default * 3.735111370627234931;
				output2[vectorPos + j / outSqu + (rowIndex % outRow / outSqu) * (outCol / outSqu)][mod(j - rowIndex, outSqu) * (output2[0].size() / outSqu) + (long long)rowIndex / outRow * outSqu + rowIndex % outSqu] = b_default * 3.735111370627234931;
				output3[vectorPos + j / outSqu + (rowIndex % outRow / outSqu) * (outCol / outSqu)][mod(j - rowIndex, outSqu) * (output3[0].size() / outSqu) + (long long)rowIndex / outRow * outSqu + rowIndex % outSqu] = b_default * 3.735111370627234931;
			}
		}
		rowIndex++;
		if (rowIndex == output1[0].size() * (outRow / outSqu) * (outCol / outSqu) / outCol * (output1[0].size() / ((long long)outSqu * outSqu) - 1) / (output1[0].size() / ((long long)outSqu * outSqu)))
		{
			rowIndex = 0;
			vectorPos += (outRow / outSqu) * (outCol / outSqu);
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < dataSize / 3; j++)
		{
			delete[] temp[i][j];
		}
		delete[] temp[i];
	}
	delete[] temp;
}

void dataCombineEncode(vector<double> output[], vector<double> input1[], vector<double> input2[], int dataSize, int rowN, int colN, int inSqu, int outSqu, vector<double> empty, int resultNum)
{
	double** temp = new double* [dataSize];
	for (int i = 0; i < dataSize; i++)
	{
		temp[i] = new double[colN] { 0 };
	}
	int inRow = ceil(rowN * 1.0 / inSqu);
	int inCol = ceil(colN * 1.0 / inSqu);
	int rowPos = slotCount * inRow * inCol * (slotCount / (inSqu * inSqu) - 1) / (slotCount / (inSqu * inSqu)) / (inCol * inSqu);
	for (int i = 0; i < dataSize / 2; i++)
	{
		for (int j = 0; j < colN; j++)
		{
			temp[i * 2][j] = input1[(i / rowPos) * inRow * inCol + (j / inSqu) * inCol + i % (inCol * inSqu) / inSqu][(long long)(i % rowPos / (inCol * inSqu)) * inSqu + i % inSqu + j % inSqu * slotCount / inSqu];
			temp[i * 2 + 1][j] = input2[(i / rowPos) * inRow * inCol + (j / inSqu) * inCol + i % (inCol * inSqu) / inSqu][(long long)(i % rowPos / (inCol * inSqu)) * inSqu + i % inSqu + j % inSqu * slotCount / inSqu];
		}
	}
	for (int i = 0; i < resultNum; i++)
	{
		output[i] = empty;
	}
	int outRow = ceil(rowN * 1.0 / outSqu) * outSqu;
	int outCol = ceil(colN * 1.0 / outSqu) * outSqu;
	int rowIndex = 0;
	int vectorPos = 0;
	for (int i = 0; i < dataSize; i++)
	{
		for (int j = 0; j < outCol; j++)
		{
			if (j < colN)
			{
				output[vectorPos + j / outSqu + (rowIndex % outRow / outSqu) * (outCol / outSqu)][mod(j - rowIndex, outSqu) * (output[0].size() / outSqu) + (long long)rowIndex / outRow * outSqu + rowIndex % outSqu] = temp[i][j];
			}
		}
		rowIndex++;
		if (rowIndex == output[0].size() * (outRow / outSqu) * (outCol / outSqu) / outCol * (output[0].size() / ((long long)outSqu * outSqu) - 1) / (output[0].size() / ((long long)outSqu * outSqu)))
		{
			rowIndex = 0;
			vectorPos += (outRow / outSqu) * (outCol / outSqu);
		}
	}
	for (int i = 0; i < dataSize; i++)
	{
		delete[] temp[i];
	}
	delete[] temp;
}

void dataCombineEncodeThree(vector<double> output[], vector<double> input1[], vector<double> input2[], vector<double> input3[], int dataSize, int rowN, int colN, int inSqu, int outSqu, vector<double> empty, int resultNum)
{
	double** temp = new double* [dataSize];
	for (int i = 0; i < dataSize; i++)
	{
		temp[i] = new double[colN] { 0 };
	}
	int inRow = ceil(rowN * 1.0 / inSqu);
	int inCol = ceil(colN * 1.0 / inSqu);
	int rowPos = slotCount * inRow * inCol * (slotCount / (inSqu * inSqu) - 1) / (slotCount / (inSqu * inSqu)) / (inCol * inSqu);
	for (int i = 0; i < dataSize / 3; i++)
	{
		for (int j = 0; j < colN; j++)
		{
			temp[i * 3][j] = input1[(i / rowPos) * inRow * inCol + (j / inSqu) * inCol + i % (inCol * inSqu) / inSqu][(long long)(i % rowPos / (inCol * inSqu)) * inSqu + i % inSqu + j % inSqu * slotCount / inSqu];
			temp[i * 3 + 1][j] = input2[(i / rowPos) * inRow * inCol + (j / inSqu) * inCol + i % (inCol * inSqu) / inSqu][(long long)(i % rowPos / (inCol * inSqu)) * inSqu + i % inSqu + j % inSqu * slotCount / inSqu];
			temp[i * 3 + 2][j] = input3[(i / rowPos) * inRow * inCol + (j / inSqu) * inCol + i % (inCol * inSqu) / inSqu][(long long)(i % rowPos / (inCol * inSqu)) * inSqu + i % inSqu + j % inSqu * slotCount / inSqu];
		}
	}
	for (int i = 0; i < resultNum; i++)
	{
		output[i] = empty;
	}
	int outRow = ceil(rowN * 1.0 / outSqu) * outSqu;
	int outCol = ceil(colN * 1.0 / outSqu) * outSqu;
	int rowIndex = 0;
	int vectorPos = 0;
	for (int i = 0; i < dataSize; i++)
	{
		for (int j = 0; j < outCol; j++)
		{
			if (j < colN)
			{
				output[vectorPos + j / outSqu + (rowIndex % outRow / outSqu) * (outCol / outSqu)][mod(j - rowIndex, outSqu) * (output[0].size() / outSqu) + (long long)rowIndex / outRow * outSqu + rowIndex % outSqu] = temp[i][j];
			}
		}
		rowIndex++;
		if (rowIndex == output[0].size() * (outRow / outSqu) * (outCol / outSqu) / outCol * (output[0].size() / ((long long)outSqu * outSqu) - 1) / (output[0].size() / ((long long)outSqu * outSqu)))
		{
			rowIndex = 0;
			vectorPos += (outRow / outSqu) * (outCol / outSqu);
		}
	}
	for (int i = 0; i < dataSize; i++)
	{
		delete[] temp[i];
	}
	delete[] temp;
}

void dataDivide(vector<double> input[], vector<double> output1[], vector<double> output2[], int dataSize, int rowN, int colN, int inSqu, int outSqu, vector<double> empty, int resultNum, double b_default)
{
	double*** temp = new double** [2];
	for (int i = 0; i < 2; i++)
	{
		temp[i] = new double* [dataSize / 2];
		for (int j = 0; j < dataSize / 2; j++)
		{
			temp[i][j] = new double[colN];
		}
	}
	int inRow = ceil(rowN * 1.0 / inSqu);
	int inCol = ceil(colN * 1.0 / inSqu);
	int rowPos = slotCount * inRow * inCol * (slotCount / (inSqu * inSqu) - 1) / (slotCount / (inSqu * inSqu)) / (inCol * inSqu);
	for (int i = 0; i < dataSize; i++)
	{
		for (int j = 0; j < colN; j++)
		{
			temp[i % 2][i / 2][j] = input[(i / rowPos) * inRow * inCol + (j / inSqu) * inCol + i % (inCol * inSqu) / inSqu][(long long)(i % rowPos / (inCol * inSqu)) * inSqu + i % inSqu + j % inSqu * slotCount / inSqu];
		}
	}
	for (int i = 0; i < resultNum; i++)
	{
		output1[i] = empty;
		output2[i] = empty;
	}
	int outRow = ceil(rowN * 1.0 / outSqu) * outSqu;
	int outCol = ceil(colN * 1.0 / outSqu) * outSqu;
	int rowIndex = 0;
	int vectorPos = 0;
	for (int i = 0; i < dataSize / 2; i++)
	{
		for (int j = 0; j < outCol; j++)
		{
			if (j < colN)
			{
				output1[vectorPos + j / outSqu * (outCol / outSqu) + (rowIndex % outRow / outSqu)][j % outSqu * (output1[0].size() / outSqu) + (long long)rowIndex / outRow * outSqu + rowIndex % outSqu] = temp[0][i][j] / pow(2, 7.75);
				output2[vectorPos + j / outSqu * (outCol / outSqu) + (rowIndex % outRow / outSqu)][j % outSqu * (output2[0].size() / outSqu) + (long long)rowIndex / outRow * outSqu + rowIndex % outSqu] = temp[1][i][j] / pow(2, 7.75);
			}
			else
			{
				output1[vectorPos + j / outSqu * (outCol / outSqu) + (rowIndex % outRow / outSqu)][j % outSqu * (output1[0].size() / outSqu) + (long long)rowIndex / outRow * outSqu + rowIndex % outSqu] = b_default * pow(2, 7.75);
				output2[vectorPos + j / outSqu * (outCol / outSqu) + (rowIndex % outRow / outSqu)][j % outSqu * (output2[0].size() / outSqu) + (long long)rowIndex / outRow * outSqu + rowIndex % outSqu] = b_default * pow(2, 7.75);
			}
		}
		rowIndex++;
		if (rowIndex == output1[0].size() * (outRow / outSqu) * (outCol / outSqu) / outCol * (output1[0].size() / ((long long)outSqu * outSqu) - 1) / (output1[0].size() / ((long long)outSqu * outSqu)))
		{
			rowIndex = 0;
			vectorPos += (outRow / outSqu) * (outCol / outSqu);
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < dataSize / 2; j++)
		{
			delete[] temp[i][j];
		}
		delete[] temp[i];
	}
	delete[] temp;
}

void dataDivideThree(vector<double> input[], vector<double> output1[], vector<double> output2[], vector<double> output3[], int dataSize, int rowN, int colN, int inSqu, int outSqu, vector<double> empty, int resultNum, double b_default)
{
	double*** temp = new double** [3];
	for (int i = 0; i < 3; i++)
	{
		temp[i] = new double* [dataSize / 3];
		for (int j = 0; j < dataSize / 3; j++)
		{
			temp[i][j] = new double[colN];
		}
	}
	int inRow = ceil(rowN * 1.0 / inSqu);
	int inCol = ceil(colN * 1.0 / inSqu);
	int rowPos = slotCount * inRow * inCol * (slotCount / (inSqu * inSqu) - 1) / (slotCount / (inSqu * inSqu)) / (inCol * inSqu);
	for (int i = 0; i < dataSize; i++)
	{
		for (int j = 0; j < colN; j++)
		{
			temp[i % 3][i / 3][j] = input[(i / rowPos) * inRow * inCol + (j / inSqu) * inCol + i % (inCol * inSqu) / inSqu][(long long)(i % rowPos / (inCol * inSqu)) * inSqu + i % inSqu + j % inSqu * slotCount / inSqu];
		}
	}
	for (int i = 0; i < resultNum; i++)
	{
		output1[i] = empty;
		output2[i] = empty;
		output3[i] = empty;
	}
	int outRow = ceil(rowN * 1.0 / outSqu) * outSqu;
	int outCol = ceil(colN * 1.0 / outSqu) * outSqu;
	int rowIndex = 0;
	int vectorPos = 0;
	for (int i = 0; i < dataSize / 3; i++)
	{
		for (int j = 0; j < outCol; j++)
		{
			if (j < colN)
			{
				output1[vectorPos + j / outSqu * (outCol / outSqu) + (rowIndex % outRow / outSqu)][j % outSqu * (output1[0].size() / outSqu) + (long long)rowIndex / outRow * outSqu + rowIndex % outSqu] = temp[0][i][j] / pow(2, 7.75);
				output2[vectorPos + j / outSqu * (outCol / outSqu) + (rowIndex % outRow / outSqu)][j % outSqu * (output2[0].size() / outSqu) + (long long)rowIndex / outRow * outSqu + rowIndex % outSqu] = temp[1][i][j] / pow(2, 7.75);
				output3[vectorPos + j / outSqu * (outCol / outSqu) + (rowIndex % outRow / outSqu)][j % outSqu * (output3[0].size() / outSqu) + (long long)rowIndex / outRow * outSqu + rowIndex % outSqu] = temp[2][i][j] / pow(2, 7.75);
			}
			else
			{
				output1[vectorPos + j / outSqu * (outCol / outSqu) + (rowIndex % outRow / outSqu)][j % outSqu * (output1[0].size() / outSqu) + (long long)rowIndex / outRow * outSqu + rowIndex % outSqu] = b_default * pow(2, 7.75);
				output2[vectorPos + j / outSqu * (outCol / outSqu) + (rowIndex % outRow / outSqu)][j % outSqu * (output2[0].size() / outSqu) + (long long)rowIndex / outRow * outSqu + rowIndex % outSqu] = b_default * pow(2, 7.75);
				output3[vectorPos + j / outSqu * (outCol / outSqu) + (rowIndex % outRow / outSqu)][j % outSqu * (output3[0].size() / outSqu) + (long long)rowIndex / outRow * outSqu + rowIndex % outSqu] = b_default * pow(2, 7.75);
			}
		}
		rowIndex++;
		if (rowIndex == output1[0].size() * (outRow / outSqu) * (outCol / outSqu) / outCol * (output1[0].size() / ((long long)outSqu * outSqu) - 1) / (output1[0].size() / ((long long)outSqu * outSqu)))
		{
			rowIndex = 0;
			vectorPos += (outRow / outSqu) * (outCol / outSqu);
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < dataSize / 3; j++)
		{
			delete[] temp[i][j];
		}
		delete[] temp[i];
	}
	delete[] temp;
}

/*测试用*/
void matrixTest(vector<double> v[], int rowN, int colN, int squN, int resultNum, int sampleCount)
{
	int finRow = ceil(rowN * 1.0 / squN);
	int finCol = ceil(colN * 1.0 / squN);
	int inputLen = resultNum / (finRow * finCol);
	int rowPos = slotCount * finRow * finCol * (slotCount / (squN * squN) - 1) / (slotCount / (squN * squN)) / (finCol * squN);
	int index1, index2, index3, index4, index5, index6;
	ofstream out("D:\\1.txt");
	for (int i = 0; i < sampleCount; i++)
	{
		index1 = i / rowPos;
		index2 = i % (finCol * squN) / squN;
		index3 = i % rowPos / (finCol * squN);
		index4 = i % squN;
		for (int j = 0; j < colN; j++)
		{
			index5 = j / squN;
			index6 = j % squN * slotCount / squN;
			out << v[index1 * finRow * finCol + index5 * finCol + index2][(long long)index3 * squN + index4 + index6] << endl;
		}
	}
	out.close();
}

int main()
{
	int updatestep[64]{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 28, 30, 31, 32, 33, 34, 36, 40, 48, 56, 60, 62, 63, 64, 65, 66, 68, 72, 80, 96, 112, 120, 124, 126, 127, 128, 129, 130, 132, 136, 144, 160, 192, 224, 240, 248, 252, 254, 255 };
	int updatezero[96]{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 46, 47, 48, 49, 50, 52, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 76, 78, 79, 80, 81, 96, 112, 120, 124, 126, 127, 128, 129, 130, 132, 136, 144, 160, 192, 224, 240, 248, 252, 254, 255 };
	//数据提供者：数据准备
	ifstream in1("D:\\TrainData.csv");
	static double data[row][col][time];
	string line;
	int a1, d1 = 0;
	getline(in1, line);
	while (getline(in1, line))
	{
		for (int i = 0; i < time; i++)
		{
			for (int j = 0; j < col; j++)
			{
				a1 = line.find(',');
				if (line.substr(0, a1) != "")
				{
					data[d1][j][i] = atof(line.substr(0, a1).c_str());
				}
				else
				{
					data[d1][j][i] = 1000000;
				}
				line.erase(0, (long long)a1 + 1);
			}
		}
		d1++;
	}
	in1.close();
	static double average[col];
	static double avect[col];
	for (int i = 0; i < col; i++)
	{
		for (int j = 0; j < row; j++)
		{
			for (int k = 0; k < time; k++)
			{
				if (data[j][i][k] != 1000000)
				{
					avect[i]++;
					average[i] += data[j][i][k];
				}
			}
		}
		average[i] = average[i] / avect[i];
	}
	static double stdev[col];
	for (int i = 0; i < col; i++)
	{
		for (int j = 0; j < row; j++)
		{
			for (int k = 0; k < time; k++)
			{
				if (data[j][i][k] != 1000000)
				{
					stdev[i] += pow(data[j][i][k] - average[i], 2);
				}
			}
		}
		stdev[i] = sqrt(stdev[i] / avect[i]);
	}
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			for (int k = 0; k < time; k++)
			{
				if (data[i][j][k] != 1000000)
				{
					data[i][j][k] = (data[i][j][k] - average[j]) / stdev[j];
				}
				else
				{
					data[i][j][k] = 0;
				}
			}
		}
	}
	ifstream in2("D:\\gsutil\\TestData\\physionet.org\\files\\challenge-2012\\1.0.0\\Outcome.csv");
	static double outcome[row];
	getline(in2, line);
	for (int i = 0; i < row; i++)
	{
		getline(in2, line);
		outcome[i] = atof(line.c_str());
	}
	in2.close();
	//研究者：训练准备
	int epoch = 20;
	static double b_0[hidenode], b_1[hidenode], b_2[hidenode], b_3[hidenode], b_4[hidenode];
	static double w_0[col][hidenode], w_1[hidenode][hidenode], w_2[hidenode][hidenode], w_3[hidenode][hidenode], w_4[hidenode][hidenode];
	static double wh_0[hidenode][hidenode], wh_1[hidenode][hidenode], wh_2[hidenode][hidenode], wh_3[hidenode][hidenode], wh_4[hidenode][hidenode];
	double b1 = 0;
	double b1_delta = 0;
	static double w1[hidenode];
	static double w1_delta[hidenode];
	static double wh_0_delta[hidenode][hidenode], wh_1_delta[hidenode][hidenode], wh_2_delta[hidenode][hidenode], wh_3_delta[hidenode][hidenode], wh_4_delta[hidenode][hidenode];
	static double w_0_delta[col][hidenode], w_1_delta[hidenode][hidenode], w_2_delta[hidenode][hidenode], w_3_delta[hidenode][hidenode], w_4_delta[hidenode][hidenode];
	static double b_0_delta[hidenode], b_1_delta[hidenode], b_2_delta[hidenode], b_3_delta[hidenode], b_4_delta[hidenode];
	double alpha;
	ifstream in3("D:\\param.bin", ios::binary);
	for (int i = 0; i < col; i++)
	{
		for (int j = 0; j < hidenode; j++)
		{
			in3.read((char*)&w_0[i][j], sizeof(double));
		}
	}
	for (int i = 0; i < hidenode; i++)
	{
		for (int j = 0; j < hidenode; j++)
		{
			in3.read((char*)&w_1[i][j], sizeof(double));
			in3.read((char*)&w_2[i][j], sizeof(double));
			in3.read((char*)&w_3[i][j], sizeof(double));
			in3.read((char*)&w_4[i][j], sizeof(double));
			in3.read((char*)&wh_0[i][j], sizeof(double));
			in3.read((char*)&wh_1[i][j], sizeof(double));
			in3.read((char*)&wh_2[i][j], sizeof(double));
			in3.read((char*)&wh_3[i][j], sizeof(double));
			in3.read((char*)&wh_4[i][j], sizeof(double));
		}
	}
	for (int i = 0; i < hidenode; i++)
	{
		in3.read((char*)&w1[i], sizeof(double));
	}
	in3.close();
	double xieLvKu[32]{ 0, 0, 0, 0, 0, 0, 7.7996707283797329e-05, 0.00021195555448024638, 0.00057570106393416728, 0.0015615719622341291, 0.0042202277676500807, 0.011293359037806703, 0.029439663215475222, 0.071777048844550773, 0.14973849934787756, 0.23105857863000490, 0.23105857863000490, 0.14973849934787742, 0.071777048844551050, 0.029439663215475087, 0.011293359037806816, 0.0042202277676500755, 0.0015615719622340540, 0.00057570106393423082, 0.00021195555448005887, 7.7996707283922717e-05, 0, 0, 0, 0, 0, 0 };
	double jieJuKu[32]{ 4.5397868702434395e-05, 4.5397868702434395e-05, 4.5397868702434395e-05, 4.5397868702434395e-05, 4.5397868702434395e-05, 4.5397868702434395e-05, 0.00082536494154040769, 0.0020309945663084493, 0.0049409586419398160, 0.011842054930039550, 0.027793989762535259, 0.063159646113318368, 0.13574486282399245, 0.26275701971121912, 0.41867992071787263, 0.5, 0.5, 0.58132007928212748, 0.73724298028878021, 0.86425513717600810, 0.93684035388668119, 0.97220601023746489, 0.98815794506996102, 0.99505904135805978, 0.99796900543369316, 0.99917463505845838, 0.99995460213129761, 0.99995460213129761, 0.99995460213129761, 0.99995460213129761, 0.99995460213129761, 0.99995460213129761 };
	//研究者：加密环境配置
	EncryptionParameters parms(scheme_type::ckks);
	parms.set_poly_modulus_degree(slotCount * 2);
	parms.set_coeff_modulus(CoeffModulus::Create(slotCount * 2, { 31, 31, 31, 31, 31, 31, 31 }));
	SEALContext context(parms);
	KeyGenerator keygen(context, partyNum);
	static SecretKey secret_key[partyNum];
	secret_key[0] = keygen.secret_key();
	for (int i = 1; i < partyNum; i++)
	{
		KeyGenerator keygenNew(context, 1);
		secret_key[i] = keygenNew.secret_key();
		keygen.skadd(secret_key[i]);
	}
	PublicKey public_key;
	keygen.create_public_key(public_key);
	RelinKeys relin_keys;
	keygen.create_relin_keys(relin_keys);
	GaloisKeys galois_keys;
	keygen.create_galois_keys(galois_keys);
	Encryptor encryptor(context, public_key);
	Evaluator evaluator(context);
	vector<Decryptor*> decryptor;
	for (int i = 0; i < partyNum; i++)
	{
		decryptor.push_back(new Decryptor(context, secret_key[i]));
	}
	CKKSEncoder ckks_encoder(context);
	/*中间变量声明*/
	int Pos;
	vector<double> empty, temp;
	Plaintext yicixiang_2Plain, yicixiang_3Plain, one_1Plain, one_2Plain, one_3Plain, one_filterPlain, bigOne_2Plain, half_3Plain, tempPlain;
	Ciphertext wb_0Cip, wb_1Cip, wb_2Cip, wb_3Cip, wb_4Cip, wh_0Cip, wh_1Cip, wh_2Cip, wh_3Cip, wh_4Cip, wb1_Cip, wb_1TranCip, wb_2TranCip, wb_3TranCip, wb_4TranCip, wh_0TranCip, wh_1TranCip, wh_2TranCip, wh_3TranCip, wh_4TranCip, wh_4_deltaCip, wh_3_deltaCip, wh_2_deltaCip, wh_1_deltaCip, wh_0_deltaCip, wb_4_deltaCip, wb_3_deltaCip, wb_2_deltaCip, wb_1_deltaCip, wb_0_deltaCip, zeroCip, tempCip;
	static vector<double> layer_0_0_0[78 * (int)1.000], layer_0_0_1[78 * (int)1.000], layer_2_error[2 * (int)1.000], xieLv[2 * (int)1.000], jieJu[2 * (int)1.000], y[2 * (int)1.000], coeff[26 * (int)1.000], noises[3][56 * (int)1.000], backnoises[3][partyNum - 1][56 * (int)1.000], temps[56 * (int)1.000];
	static Plaintext wb_0filterPlain[6], whb_filterPlain[4], wb1_filterPlain[2], layer_2_errorPlain[2 * (int)1.000], xieLvPlain[2 * (int)1.000], jieJuPlain[2 * (int)1.000], pDecPlain[52 * (int)1.000][partyNum - 1], tempsPlain[26 * (int)1.000];
	static Ciphertext wb_0EncodeCip[6], wb_1EncodeCip[4], wb_2EncodeCip[4], wb_3EncodeCip[4], wb_4EncodeCip[4], wh_0EncodeCip[4], wh_1EncodeCip[4], wh_2EncodeCip[4], wh_3EncodeCip[4], wh_4EncodeCip[4], wb1_EncodeCip[2], wb1_NoiseCip[partyNum - 1], wh_4_NoiseCip[partyNum - 1], wh_4_TranNoiseCip[partyNum - 1], wh_3_NoiseCip[partyNum - 1], wh_3_TranNoiseCip[partyNum - 1], wh_2_NoiseCip[partyNum - 1], wh_2_TranNoiseCip[partyNum - 1], wh_1_NoiseCip[partyNum - 1], wh_1_TranNoiseCip[partyNum - 1], wh_0_NoiseCip[partyNum - 1], wh_0_TranNoiseCip[partyNum - 1], wb_4_NoiseCip[partyNum - 1], wb_4_TranNoiseCip[partyNum - 1], wb_3_NoiseCip[partyNum - 1], wb_3_TranNoiseCip[partyNum - 1], wb_2_NoiseCip[partyNum - 1], wb_2_TranNoiseCip[partyNum - 1], wb_1_NoiseCip[partyNum - 1], wb_1_TranNoiseCip[partyNum - 1], wb_0_NoiseCip[partyNum - 1], o1_0_0Cip[52 * (int)1.000], o1_0_1Cip[52 * (int)1.000], o1_1_0Cip[28 * (int)1.000], o1_1_1Cip[28 * (int)1.000], o1_2_0Cip[16 * (int)1.000], o1_2_1Cip[16 * (int)1.000], o1_3_0Cip[8 * (int)1.000], o1_3_1Cip[8 * (int)1.000], o1_4_0Cip[4 * (int)1.000], o1_4_1Cip[4 * (int)1.000], o1_4_2Cip[4 * (int)1.000], layer_0_0Cip[78 * (int)1.000], layer_1_0_0Cip[52 * (int)1.000], layer_1_0_0OriCip[52 * (int)1.000], layer_1_0_1Cip[56 * (int)1.000], layer_1_0_1OriCip[56 * (int)1.000], layer_1_1_0Cip[28 * (int)1.000], layer_1_1_0OriCip[28 * (int)1.000], layer_1_1_1Cip[32 * (int)1.000], layer_1_1_1OriCip[32 * (int)1.000], layer_1_2_0Cip[16 * (int)1.000], layer_1_2_0OriCip[16 * (int)1.000], layer_1_2_1Cip[16 * (int)1.000], layer_1_2_1OriCip[16 * (int)1.000], layer_1_3_0Cip[8 * (int)1.000], layer_1_3_0OriCip[8 * (int)1.000], layer_1_3_1Cip[12 * (int)1.000], layer_1_3_1OriCip[12 * (int)1.000], layer_1_4_0Cip[4 * (int)1.000], layer_1_4_0OriCip[4 * (int)1.000], layer_1_4_1Cip[4 * (int)1.000], layer_1_4_1OriCip[4 * (int)1.000], layer_1_4_2Cip[4 * (int)1.000], layer_2Cip[2 * (int)1.000], layer_2_deltaCip[2 * (int)1.000], layer_1_4_deltaCip[3][4 * (int)1.000], layer_1_4_deltaOriCip[3][4 * (int)1.000][16], layer_1_3_deltaCip[2][8 * (int)1.000], layer_1_3_deltaOriCip[2][8 * (int)1.000][16], layer_1_2_deltaCip[2][16 * (int)1.000], layer_1_2_deltaOriCip[2][16 * (int)1.000][16], layer_1_1_deltaCip[2][28 * (int)1.000], layer_1_1_deltaOriCip[2][28 * (int)1.000][16], layer_1_0_deltaCip[2][52 * (int)1.000], layer_1_0_deltaOriCip[2][52 * (int)1.000][16], paramCip[2][96], noisesCip[partyNum - 1][56 * (int)1.000], noisesOriCip[partyNum - 1][56 * (int)1.000], tempsCip[4992 * (int)1.000];
	/*初始化参数密文*/
	ckks_encoder.encode(0, levelZeroScale, tempPlain);
	encryptor.encrypt(tempPlain, zeroCip);
	wb_0Cip = zeroCip;
	wb_1Cip = zeroCip;
	wb_2Cip = zeroCip;
	wb_3Cip = zeroCip;
	wb_4Cip = zeroCip;
	wh_0Cip = zeroCip;
	wh_1Cip = zeroCip;
	wh_2Cip = zeroCip;
	wh_3Cip = zeroCip;
	wh_4Cip = zeroCip;
	wb1_Cip = zeroCip;
	wb_1TranCip = zeroCip;
	wb_2TranCip = zeroCip;
	wb_3TranCip = zeroCip;
	wb_4TranCip = zeroCip;
	wh_0TranCip = zeroCip;
	wh_1TranCip = zeroCip;
	wh_2TranCip = zeroCip;
	wh_3TranCip = zeroCip;
	wh_4TranCip = zeroCip;
	for (int i = 1; i < partyNum; i++)
	{
		wb1_NoiseCip[i - 1] = zeroCip;
		wh_4_NoiseCip[i - 1] = zeroCip;
		wh_4_TranNoiseCip[i - 1] = zeroCip;
		wh_3_NoiseCip[i - 1] = zeroCip;
		wh_3_TranNoiseCip[i - 1] = zeroCip;
		wh_2_NoiseCip[i - 1] = zeroCip;
		wh_2_TranNoiseCip[i - 1] = zeroCip;
		wh_1_NoiseCip[i - 1] = zeroCip;
		wh_1_TranNoiseCip[i - 1] = zeroCip;
		wh_0_NoiseCip[i - 1] = zeroCip;
		wh_0_TranNoiseCip[i - 1] = zeroCip;
		wb_4_NoiseCip[i - 1] = zeroCip;
		wb_4_TranNoiseCip[i - 1] = zeroCip;
		wb_3_NoiseCip[i - 1] = zeroCip;
		wb_3_TranNoiseCip[i - 1] = zeroCip;
		wb_2_NoiseCip[i - 1] = zeroCip;
		wb_2_TranNoiseCip[i - 1] = zeroCip;
		wb_1_NoiseCip[i - 1] = zeroCip;
		wb_1_TranNoiseCip[i - 1] = zeroCip;
		wb_0_NoiseCip[i - 1] = zeroCip;
	}
	/*设定空白向量*/
	for (int i = 0; i < slotCount; i++)
	{
		empty.push_back(0);
	}
	for (int i = 0; i < 56 * (int)1.000; i++)
	{
		temps[i] = empty;
		for (int j = 0; j < 3; j++)
		{
			noises[j][i] = empty;
			for (int k = 0; k < partyNum - 1; k++)
			{
				backnoises[j][k][i] = empty;
			}
		}
	}
	/*设定参数filter及其他杂项*/
	for (int i = 0; i < 6; i++)
	{
		temp = empty;
		for (int j = 0; j < pow(2, ceil(log2(hidenode * col / 6) / 2)); j++)
		{
			for (int k = 0; k < pow(2, ceil(log2(hidenode * col / 6) / 2)) * 2; k++)
			{
				temp[(long long)j * slotCount / pow(2, ceil(log2(hidenode * col / 6) / 2)) + k + i * (pow(2, ceil(log2(hidenode * col / 6) / 2)) * 2)] = 1;
			}
		}
		ckks_encoder.encode(temp, levelZeroScale, wb_0filterPlain[i]);
	}
	for (int i = 0; i < 4; i++)
	{
		temp = empty;
		for (int j = 0; j < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); j++)
		{
			for (int k = 0; k < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4; k++)
			{
				temp[(long long)j * slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + i * (pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4)] = 1;
			}
		}
		ckks_encoder.encode(temp, levelZeroScale, whb_filterPlain[i]);
	}
	for (int i = 0; i < 2 * (int)1.000; i++)
	{
		layer_2_error[i] = empty;
	}
	for (int i = 0; i < 2; i++)
	{
		temp = empty;
		for (int j = 0; j < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); j++)
		{
			for (int k = 0; k < (slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) - pow(2, ceil(log2(hidenode * hidenode / 4) / 2))) / 2; k++)
			{
				temp[(long long)j * slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + (long long)i * slotCount / pow(2, 1 + ceil(log2(hidenode * hidenode / 4) / 2))] = 1;
			}
		}
		ckks_encoder.encode(temp, levelZeroScale, wb1_filterPlain[i]);
	}
	temp = empty;
	for (int i = 0; i < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); i++)
	{
		temp[(long long)i * slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2))] = 1;
	}
	ckks_encoder.encode(temp, levelThree, levelThreeScale, one_filterPlain);
	temp = empty;
	for (int i = (pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) - 1) * slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); i < (pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) - 1) * slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) - pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); i++)
	{
		temp[i] = pow(2, 15.5);
	}
	ckks_encoder.encode(temp, levelThree, levelThreeScale, one_3Plain);
	ckks_encoder.encode(1, levelOne, levelOneScale, one_1Plain);
	ckks_encoder.encode(1, levelTwo, levelTwoScale, one_2Plain);
	ckks_encoder.encode(pow(2, 15.5 / 3.0), levelTwo, levelTwoScale, bigOne_2Plain);
	ckks_encoder.encode(1660.840113212959707599, levelTwo, levelTwoScale, yicixiang_2Plain);
	ckks_encoder.encode(1660.840113212959707599, levelThree, levelThreeScale, yicixiang_3Plain);
	ckks_encoder.encode(830.420056606479856743, levelThree, levelThreeScale, half_3Plain);
	/*开始训练*/
	for (int i = 0; i < epoch * row / batchsize; i++)
	{
		alpha = 0.001 * pow(0.5, (i / (row / batchsize)) / 100.0);
		//研究者：参数加密
		if (i == 0)
		{
			parameterEncode(wb_0Cip, hidenode, col, pow(2, ceil(log2(hidenode * col / 6) / 2)), evaluator, ckks_encoder, empty, w_0, b_0);
			parameterEncode(wb_1Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w_1, b_1);
			parameterEncode(wb_2Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w_2, b_2);
			parameterEncode(wb_3Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w_3, b_3);
			parameterEncode(wb_4Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w_4, b_4);
			parameterEncode(wh_0Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_0);
			parameterEncode(wh_1Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_1);
			parameterEncode(wh_2Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_2);
			parameterEncode(wh_3Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_3);
			parameterEncode(wh_4Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_4);
			outParameterEncode(wb1_Cip, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w1, b1);
			parameterEncodeTranspose(wh_4TranCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_4);
			parameterEncodeTranspose(wh_3TranCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_3);
			parameterEncodeTranspose(wh_2TranCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_2);
			parameterEncodeTranspose(wh_1TranCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_1);
			parameterEncodeTranspose(wh_0TranCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_0);
			parameterEncodeTranspose(wb_4TranCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w_4);
			parameterEncodeTranspose(wb_3TranCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w_3);
			parameterEncodeTranspose(wb_2TranCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w_2);
			parameterEncodeTranspose(wb_1TranCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w_1);
		}
		//研究者：将10个密文参数发送
		//数据提供者：正向传播1,3,5,7...,47
		for (int j = 0; j < 6; j++)
		{
			evaluator.multiply_plain(wb_0Cip, wb_0filterPlain[j], wb_0EncodeCip[j]);
			for (int k = pow(2, ceil(log2(hidenode * col / 6) / 2)) * 2; k < slotCount / pow(2, ceil(log2(hidenode * col / 6) / 2)); k *= 2)
			{
				evaluator.rotate_vector(wb_0EncodeCip[j], -k, galois_keys, tempCip);
				evaluator.add_inplace(wb_0EncodeCip[j], tempCip);
			}
			if (j != 0)
			{
				evaluator.rotate_vector_inplace(wb_0EncodeCip[j], j * pow(2, ceil(log2(hidenode * col / 6) / 2)) * 2, galois_keys);
			}
			evaluator.rescale_to_next_inplace(wb_0EncodeCip[j]);
		}
		parameterRotate(wb_0EncodeCip, hidenode, col, pow(2, ceil(log2(hidenode * col / 6) / 2)), paramCip[0], evaluator, galois_keys, one_1Plain);
		for (int j = 0; j < 78 * (int)1.000; j++)
		{
			layer_0_0_0[j] = empty;
		}
		dataEncode(layer_0_0_0, hidenode, col, pow(2, ceil(log2(hidenode * col / 6) / 2)), data, 0, i);
		matrixMultiple_CipPlain(paramCip[0], layer_0_0_0, hidenode, col, pow(2, ceil(log2(hidenode * col / 6) / 2)), 52 * (int)1.000, o1_0_0Cip, evaluator, ckks_encoder);
		for (int j = 0; j < 52 * (int)1.000; j++)
		{
			evaluator.rescale_to_next_inplace(o1_0_0Cip[j]);
			evaluator.square(o1_0_0Cip[j], layer_1_0_0Cip[j]);
			evaluator.relinearize_inplace(layer_1_0_0Cip[j], relin_keys);
			evaluator.multiply_plain(o1_0_0Cip[j], yicixiang_3Plain, tempCip);
			evaluator.add_inplace(layer_1_0_0Cip[j], tempCip);
			evaluator.rescale_to_next_inplace(layer_1_0_0Cip[j]);
		}
		//数据提供者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < 52 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					noises[j - 1][k][m] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
				}
				ckks_encoder.encode(noises[j - 1][k], levelFour, levelFourScale, tempPlain);
				decryptor[j]->partial_decrypt(layer_1_0_0Cip[k], tempPlain, pDecPlain[k][j - 1]);
				for (int m = 0; m < slotCount; m++)
				{
					noises[j - 1][k][m] /= -pow(2, 7.75);
				}
				ckks_encoder.encode(noises[j - 1][k], levelOne, levelOneScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesOriCip[j - 1][k]);
			}
			matrixRefresh(noises[j - 1], temps, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 52 * (int)1.000, batchsize * time / 2);
			for (int k = 0; k < 52 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					temps[k][m] *= 3.735111370627234931 / pow(2, 7.75);
				}
				ckks_encoder.encode(temps[k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesCip[j - 1][k]);
			}
		}
		//数据提供者：将52个密文的后半部分广播/前半部分发送，52×2×(partyNum - 1)个噪声密文发送，52×(partyNum - 1)个部分解密发送
		//研究者：完全解密并重新编码(将52个layer_1_0_0OriCip发送给数据提供者)
		for (int j = 0; j < 52 * (int)1.000; j++)
		{
			decryptor[0]->final_decrypt(layer_1_0_0Cip[j], pDecPlain[j], partyNum - 1, tempPlain);
			ckks_encoder.decode(tempPlain, noises[0][j]);
		}
		for (int j = 0; j < 52 * (int)1.000; j++)
		{
			for (int k = 0; k < slotCount; k++)
			{
				noises[0][j][k] /= pow(2, 7.75);
			}
			ckks_encoder.encode(noises[0][j], levelOne, levelOneScale, tempPlain);
			layer_1_0_0OriCip[j] = noisesOriCip[0][j];
			evaluator.add_plain_inplace(layer_1_0_0OriCip[j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_0_0OriCip[j], noisesOriCip[k - 1][j]);
				k++;
			}
		}
		matrixRefresh(noises[0], temps, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 52 * (int)1.000, batchsize * time / 2);
		for (int j = 0; j < 52 * (int)1.000; j++)
		{
			for (int k = 0; k < slotCount; k++)
			{
				temps[j][k] *= 3.735111370627234931 / pow(2, 7.75);
			}
			ckks_encoder.encode(temps[j], levelTwo, levelTwoScale, tempPlain);
			layer_1_0_0Cip[j] = noisesCip[0][j];
			evaluator.add_plain_inplace(layer_1_0_0Cip[j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_0_0Cip[j], noisesCip[k - 1][j]);
				k++;
			}
		}
		//研究者：计算wh_0 * layer_1_0_0
		for (int j = 0; j < 4; j++)
		{
			evaluator.multiply_plain(wh_0Cip, whb_filterPlain[j], wh_0EncodeCip[j]);
			for (int k = pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(wh_0EncodeCip[j], -k, galois_keys, tempCip);
				evaluator.add_inplace(wh_0EncodeCip[j], tempCip);
			}
			if (j != 0)
			{
				evaluator.rotate_vector_inplace(wh_0EncodeCip[j], j * pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4, galois_keys);
			}
			evaluator.rescale_to_next_inplace(wh_0EncodeCip[j]);
		}
		parameterRotate(wh_0EncodeCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), paramCip[1], evaluator, galois_keys, one_1Plain);
		matrixMultiple_CipCip(paramCip[1], layer_1_0_0Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 52 * (int)1.000, layer_1_0_1Cip, evaluator, galois_keys);
		for (int j = 0; j < 52 * (int)1.000; j++)
		{
			evaluator.relinearize_inplace(layer_1_0_1Cip[j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_0_1Cip[j]);
		}
		//研究者：将52个密文发送
		//数据提供者：正向传播2,4,6,8,...,48
		for (int j = 0; j < 78 * (int)1.000; j++)
		{
			layer_0_0_1[j] = empty;
		}
		dataEncode(layer_0_0_1, hidenode, col, pow(2, ceil(log2(hidenode * col / 6) / 2)), data, 1, i);
		matrixMultiple_CipPlain(paramCip[0], layer_0_0_1, hidenode, col, pow(2, ceil(log2(hidenode * col / 6) / 2)), 52 * (int)1.000, o1_0_1Cip, evaluator, ckks_encoder);
		for (int j = 0; j < 52 * (int)1.000; j++)
		{
			evaluator.rescale_to_next_inplace(o1_0_1Cip[j]);
			evaluator.add_inplace(o1_0_1Cip[j], layer_1_0_1Cip[j]);
			evaluator.square(o1_0_1Cip[j], layer_1_0_1Cip[j]);
			evaluator.relinearize_inplace(layer_1_0_1Cip[j], relin_keys);
			evaluator.multiply_plain(o1_0_1Cip[j], yicixiang_3Plain, tempCip);
			evaluator.add_inplace(layer_1_0_1Cip[j], tempCip);
			evaluator.rescale_to_next_inplace(layer_1_0_1Cip[j]);
		}
		//数据提供者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < 52 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					noises[j - 1][k][m] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
				}
				ckks_encoder.encode(noises[j - 1][k], levelFour, levelFourScale, tempPlain);
				decryptor[j]->partial_decrypt(layer_1_0_1Cip[k], tempPlain, pDecPlain[k][j - 1]);
				for (int m = 0; m < slotCount; m++)
				{
					noises[j - 1][k][m] *= -1;
				}
			}
			dataDivide(noises[j - 1], temps, temps + 28 * (int)1.000, batchsize * time / 2, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), empty, 28 * (int)1.000, 0);
			for (int k = 0; k < 56 * (int)1.000; k++)
			{
				ckks_encoder.encode(temps[k], levelOne, levelOneScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesOriCip[j - 1][k]);
			}
			dataDivideEncode(noises[j - 1], temps, temps + 28 * (int)1.000, batchsize * time / 2, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), empty, 28 * (int)1.000, 0);
			for (int k = 0; k < 56 * (int)1.000; k++)
			{
				ckks_encoder.encode(temps[k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesCip[j - 1][k]);
			}
		}
		//数据提供者：将52个密文的后半部分广播/前半部分发送，56×2×(partyNum - 1)个噪声密文发送，52×(partyNum - 1)个部分解密发送
		//研究者：完全解密并重新编码(将56个layer_1_0_1OriCip发送给数据提供者)
		for (int j = 0; j < 52 * (int)1.000; j++)
		{
			decryptor[0]->final_decrypt(layer_1_0_1Cip[j], pDecPlain[j], partyNum - 1, tempPlain);
			ckks_encoder.decode(tempPlain, noises[0][j]);
		}
		dataDivide(noises[0], temps, temps + 28 * (int)1.000, batchsize * time / 2, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), empty, 28 * (int)1.000, 1);
		for (int j = 0; j < 56 * (int)1.000; j++)
		{
			ckks_encoder.encode(temps[j], levelOne, levelOneScale, tempPlain);
			layer_1_0_1OriCip[j] = noisesOriCip[0][j];
			evaluator.add_plain_inplace(layer_1_0_1OriCip[j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_0_1OriCip[j], noisesOriCip[k - 1][j]);
				k++;
			}
		}
		dataDivideEncode(noises[0], temps, temps + 28 * (int)1.000, batchsize * time / 2, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), empty, 28 * (int)1.000, 1);
		for (int j = 0; j < 56 * (int)1.000; j++)
		{
			ckks_encoder.encode(temps[j], levelTwo, levelTwoScale, tempPlain);
			layer_1_0_1Cip[j] = noisesCip[0][j];
			evaluator.add_plain_inplace(layer_1_0_1Cip[j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_0_1Cip[j], noisesCip[k - 1][j]);
				k++;
			}
		}
		//研究者：正向传播2,6,10,14,...,46
		for (int j = 0; j < 4; j++)
		{
			evaluator.multiply_plain(wb_1Cip, whb_filterPlain[j], wb_1EncodeCip[j]);
			for (int k = pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(wb_1EncodeCip[j], -k, galois_keys, tempCip);
				evaluator.add_inplace(wb_1EncodeCip[j], tempCip);
			}
			if (j != 0)
			{
				evaluator.rotate_vector_inplace(wb_1EncodeCip[j], j * pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4, galois_keys);
			}
			evaluator.rescale_to_next_inplace(wb_1EncodeCip[j]);
		}
		parameterRotate(wb_1EncodeCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), paramCip[0], evaluator, galois_keys, one_1Plain);
		matrixMultiple_CipCip(paramCip[0], layer_1_0_1Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 28 * (int)1.000, o1_1_0Cip, evaluator, galois_keys);
		for (int j = 0; j < 28 * (int)1.000; j++)
		{
			evaluator.relinearize_inplace(o1_1_0Cip[j], relin_keys);
			evaluator.rescale_to_next_inplace(o1_1_0Cip[j]);
			evaluator.square(o1_1_0Cip[j], layer_1_1_0Cip[j]);
			evaluator.relinearize_inplace(layer_1_1_0Cip[j], relin_keys);
			evaluator.multiply_plain(o1_1_0Cip[j], yicixiang_3Plain, tempCip);
			evaluator.add_inplace(layer_1_1_0Cip[j], tempCip);
			evaluator.rescale_to_next_inplace(layer_1_1_0Cip[j]);
		}
		//研究者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < 28 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					noises[j - 1][k][m] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
				}
				ckks_encoder.encode(noises[j - 1][k], levelFour, levelFourScale, tempPlain);
				decryptor[j]->partial_decrypt(layer_1_1_0Cip[k], tempPlain, pDecPlain[k][j - 1]);
				for (int m = 0; m < slotCount; m++)
				{
					noises[j - 1][k][m] /= -pow(2, 7.75);
				}
				ckks_encoder.encode(noises[j - 1][k], levelOne, levelOneScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesOriCip[j - 1][k]);
			}
			matrixRefresh(noises[j - 1], temps, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 28 * (int)1.000, batchsize * time / 4);
			for (int k = 0; k < 28 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					temps[k][m] *= 3.735111370627234931 / pow(2, 7.75);
				}
				ckks_encoder.encode(temps[k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesCip[j - 1][k]);
			}
		}
		//研究者：将28个密文的后半部分广播/前半部分发送，28×2×(partyNum - 1)个噪声密文发送，28×(partyNum - 1)个部分解密发送
		//数据提供者：完全解密并重新编码
		for (int j = 0; j < 28 * (int)1.000; j++)
		{
			decryptor[0]->final_decrypt(layer_1_1_0Cip[j], pDecPlain[j], partyNum - 1, tempPlain);
			ckks_encoder.decode(tempPlain, noises[0][j]);
		}
		for (int j = 0; j < 28 * (int)1.000; j++)
		{
			for (int k = 0; k < slotCount; k++)
			{
				noises[0][j][k] /= pow(2, 7.75);
			}
			ckks_encoder.encode(noises[0][j], levelOne, levelOneScale, tempPlain);
			layer_1_1_0OriCip[j] = noisesOriCip[0][j];
			evaluator.add_plain_inplace(layer_1_1_0OriCip[j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_1_0OriCip[j], noisesOriCip[k - 1][j]);
				k++;
			}
		}
		matrixRefresh(noises[0], temps, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 28 * (int)1.000, batchsize * time / 4);
		for (int j = 0; j < 28 * (int)1.000; j++)
		{
			for (int k = 0; k < slotCount; k++)
			{
				temps[j][k] *= 3.735111370627234931 / pow(2, 7.75);
			}
			ckks_encoder.encode(temps[j], levelTwo, levelTwoScale, tempPlain);
			layer_1_1_0Cip[j] = noisesCip[0][j];
			evaluator.add_plain_inplace(layer_1_1_0Cip[j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_1_0Cip[j], noisesCip[k - 1][j]);
				k++;
			}
		}
		//数据提供者：计算wh_1 * layer_1_1_0
		for (int j = 0; j < 4; j++)
		{
			evaluator.multiply_plain(wh_1Cip, whb_filterPlain[j], wh_1EncodeCip[j]);
			for (int k = pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(wh_1EncodeCip[j], -k, galois_keys, tempCip);
				evaluator.add_inplace(wh_1EncodeCip[j], tempCip);
			}
			if (j != 0)
			{
				evaluator.rotate_vector_inplace(wh_1EncodeCip[j], j * pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4, galois_keys);
			}
			evaluator.rescale_to_next_inplace(wh_1EncodeCip[j]);
		}
		parameterRotate(wh_1EncodeCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), paramCip[1], evaluator, galois_keys, one_1Plain);
		matrixMultiple_CipCip(paramCip[1], layer_1_1_0Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 28 * (int)1.000, layer_1_1_1Cip, evaluator, galois_keys);
		for (int j = 0; j < 28 * (int)1.000; j++)
		{
			evaluator.relinearize_inplace(layer_1_1_1Cip[j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_1_1Cip[j]);
		}
		//数据提供者：将28个密文发送
		//研究者：正向传播4,8,12,16,...,48
		matrixMultiple_CipCip(paramCip[0], layer_1_0_1Cip + 28 * (int)1.000, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 28 * (int)1.000, o1_1_1Cip, evaluator, galois_keys);
		for (int j = 0; j < 28 * (int)1.000; j++)
		{
			evaluator.relinearize_inplace(o1_1_1Cip[j], relin_keys);
			evaluator.rescale_to_next_inplace(o1_1_1Cip[j]);
			evaluator.add_inplace(o1_1_1Cip[j], layer_1_1_1Cip[j]);
			evaluator.square(o1_1_1Cip[j], layer_1_1_1Cip[j]);
			evaluator.relinearize_inplace(layer_1_1_1Cip[j], relin_keys);
			evaluator.multiply_plain(o1_1_1Cip[j], yicixiang_3Plain, tempCip);
			evaluator.add_inplace(layer_1_1_1Cip[j], tempCip);
			evaluator.rescale_to_next_inplace(layer_1_1_1Cip[j]);
		}
		//研究者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < 28 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					noises[j - 1][k][m] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
				}
				ckks_encoder.encode(noises[j - 1][k], levelFour, levelFourScale, tempPlain);
				decryptor[j]->partial_decrypt(layer_1_1_1Cip[k], tempPlain, pDecPlain[k][j - 1]);
				for (int m = 0; m < slotCount; m++)
				{
					noises[j - 1][k][m] *= -1;
				}
			}
			dataDivide(noises[j - 1], temps, temps + 16 * (int)1.000, batchsize * time / 4, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), empty, 16 * (int)1.000, 0);
			for (int k = 0; k < 32 * (int)1.000; k++)
			{
				ckks_encoder.encode(temps[k], levelOne, levelOneScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesOriCip[j - 1][k]);
			}
			dataDivideEncode(noises[j - 1], temps, temps + 16 * (int)1.000, batchsize * time / 4, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), empty, 16 * (int)1.000, 0);
			for (int k = 0; k < 32 * (int)1.000; k++)
			{
				ckks_encoder.encode(temps[k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesCip[j - 1][k]);
			}
		}
		//研究者：将28个密文的后半部分广播/前半部分发送，32×2×(partyNum - 1)个噪声密文发送，28×(partyNum - 1)个部分解密发送
		//数据提供者：完全解密并重新编码
		for (int j = 0; j < 28 * (int)1.000; j++)
		{
			decryptor[0]->final_decrypt(layer_1_1_1Cip[j], pDecPlain[j], partyNum - 1, tempPlain);
			ckks_encoder.decode(tempPlain, noises[0][j]);
		}
		dataDivide(noises[0], temps, temps + 16 * (int)1.000, batchsize * time / 4, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), empty, 16 * (int)1.000, 1);
		for (int j = 0; j < 32 * (int)1.000; j++)
		{
			ckks_encoder.encode(temps[j], levelOne, levelOneScale, tempPlain);
			layer_1_1_1OriCip[j] = noisesOriCip[0][j];
			evaluator.add_plain_inplace(layer_1_1_1OriCip[j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_1_1OriCip[j], noisesOriCip[k - 1][j]);
				k++;
			}
		}
		dataDivideEncode(noises[0], temps, temps + 16 * (int)1.000, batchsize * time / 4, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), empty, 16 * (int)1.000, 1);
		for (int j = 0; j < 32 * (int)1.000; j++)
		{
			ckks_encoder.encode(temps[j], levelTwo, levelTwoScale, tempPlain);
			layer_1_1_1Cip[j] = noisesCip[0][j];
			evaluator.add_plain_inplace(layer_1_1_1Cip[j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_1_1Cip[j], noisesCip[k - 1][j]);
				k++;
			}
		}
		//数据提供者：正向传播4,12,20,28,36,44
		for (int j = 0; j < 4; j++)
		{
			evaluator.multiply_plain(wb_2Cip, whb_filterPlain[j], wb_2EncodeCip[j]);
			for (int k = pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(wb_2EncodeCip[j], -k, galois_keys, tempCip);
				evaluator.add_inplace(wb_2EncodeCip[j], tempCip);
			}
			if (j != 0)
			{
				evaluator.rotate_vector_inplace(wb_2EncodeCip[j], j * pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4, galois_keys);
			}
			evaluator.rescale_to_next_inplace(wb_2EncodeCip[j]);
		}
		parameterRotate(wb_2EncodeCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), paramCip[0], evaluator, galois_keys, one_1Plain);
		matrixMultiple_CipCip(paramCip[0], layer_1_1_1Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 16 * (int)1.000, o1_2_0Cip, evaluator, galois_keys);
		for (int j = 0; j < 16 * (int)1.000; j++)
		{
			evaluator.relinearize_inplace(o1_2_0Cip[j], relin_keys);
			evaluator.rescale_to_next_inplace(o1_2_0Cip[j]);
			evaluator.square(o1_2_0Cip[j], layer_1_2_0Cip[j]);
			evaluator.relinearize_inplace(layer_1_2_0Cip[j], relin_keys);
			evaluator.multiply_plain(o1_2_0Cip[j], yicixiang_3Plain, tempCip);
			evaluator.add_inplace(layer_1_2_0Cip[j], tempCip);
			evaluator.rescale_to_next_inplace(layer_1_2_0Cip[j]);
		}
		//数据提供者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < 16 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					noises[j - 1][k][m] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
				}
				ckks_encoder.encode(noises[j - 1][k], levelFour, levelFourScale, tempPlain);
				decryptor[j]->partial_decrypt(layer_1_2_0Cip[k], tempPlain, pDecPlain[k][j - 1]);
				for (int m = 0; m < slotCount; m++)
				{
					noises[j - 1][k][m] /= -pow(2, 7.75);
				}
				ckks_encoder.encode(noises[j - 1][k], levelOne, levelOneScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesOriCip[j - 1][k]);
			}
			matrixRefresh(noises[j - 1], temps, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 16 * (int)1.000, batchsize * time / 8);
			for (int k = 0; k < 16 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					temps[k][m] *= 3.735111370627234931 / pow(2, 7.75);
				}
				ckks_encoder.encode(temps[k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesCip[j - 1][k]);
			}
		}
		//数据提供者：将16个密文的后半部分广播/前半部分发送，16×2×(partyNum - 1)个噪声密文发送，16×(partyNum - 1)个部分解密发送
		//研究者：完全解密并重新编码(将16个layer_1_2_0OriCip发送给数据提供者)
		for (int j = 0; j < 16 * (int)1.000; j++)
		{
			decryptor[0]->final_decrypt(layer_1_2_0Cip[j], pDecPlain[j], partyNum - 1, tempPlain);
			ckks_encoder.decode(tempPlain, noises[0][j]);
		}
		for (int j = 0; j < 16 * (int)1.000; j++)
		{
			for (int k = 0; k < slotCount; k++)
			{
				noises[0][j][k] /= pow(2, 7.75);
			}
			ckks_encoder.encode(noises[0][j], levelOne, levelOneScale, tempPlain);
			layer_1_2_0OriCip[j] = noisesOriCip[0][j];
			evaluator.add_plain_inplace(layer_1_2_0OriCip[j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_2_0OriCip[j], noisesOriCip[k - 1][j]);
				k++;
			}
		}
		matrixRefresh(noises[0], temps, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 16 * (int)1.000, batchsize * time / 8);
		for (int j = 0; j < 16 * (int)1.000; j++)
		{
			for (int k = 0; k < slotCount; k++)
			{
				temps[j][k] *= 3.735111370627234931 / pow(2, 7.75);
			}
			ckks_encoder.encode(temps[j], levelTwo, levelTwoScale, tempPlain);
			layer_1_2_0Cip[j] = noisesCip[0][j];
			evaluator.add_plain_inplace(layer_1_2_0Cip[j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_2_0Cip[j], noisesCip[k - 1][j]);
				k++;
			}
		}
		//研究者：计算wh_2 * layer_1_2_0
		for (int j = 0; j < 4; j++)
		{
			evaluator.multiply_plain(wh_2Cip, whb_filterPlain[j], wh_2EncodeCip[j]);
			for (int k = pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(wh_2EncodeCip[j], -k, galois_keys, tempCip);
				evaluator.add_inplace(wh_2EncodeCip[j], tempCip);
			}
			if (j != 0)
			{
				evaluator.rotate_vector_inplace(wh_2EncodeCip[j], j * pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4, galois_keys);
			}
			evaluator.rescale_to_next_inplace(wh_2EncodeCip[j]);
		}
		parameterRotate(wh_2EncodeCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), paramCip[1], evaluator, galois_keys, one_1Plain);
		matrixMultiple_CipCip(paramCip[1], layer_1_2_0Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 16 * (int)1.000, layer_1_2_1Cip, evaluator, galois_keys);
		for (int j = 0; j < 16 * (int)1.000; j++)
		{
			evaluator.relinearize_inplace(layer_1_2_1Cip[j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_2_1Cip[j]);
		}
		//研究者：将16个密文发送
		//数据提供者：正向传播8,16,24,32,40,48
		matrixMultiple_CipCip(paramCip[0], layer_1_1_1Cip + 16 * (int)1.000, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 16 * (int)1.000, o1_2_1Cip, evaluator, galois_keys);
		for (int j = 0; j < 16 * (int)1.000; j++)
		{
			evaluator.relinearize_inplace(o1_2_1Cip[j], relin_keys);
			evaluator.rescale_to_next_inplace(o1_2_1Cip[j]);
			evaluator.add_inplace(o1_2_1Cip[j], layer_1_2_1Cip[j]);
			evaluator.square(o1_2_1Cip[j], layer_1_2_1Cip[j]);
			evaluator.relinearize_inplace(layer_1_2_1Cip[j], relin_keys);
			evaluator.multiply_plain(o1_2_1Cip[j], yicixiang_3Plain, tempCip);
			evaluator.add_inplace(layer_1_2_1Cip[j], tempCip);
			evaluator.rescale_to_next_inplace(layer_1_2_1Cip[j]);
		}
		//数据提供者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < 16 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					noises[j - 1][k][m] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
				}
				ckks_encoder.encode(noises[j - 1][k], levelFour, levelFourScale, tempPlain);
				decryptor[j]->partial_decrypt(layer_1_2_1Cip[k], tempPlain, pDecPlain[k][j - 1]);
				for (int m = 0; m < slotCount; m++)
				{
					noises[j - 1][k][m] *= -1;
				}
			}
			dataDivide(noises[j - 1], temps, temps + 8 * (int)1.000, batchsize * time / 8, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), empty, 8 * (int)1.000, 0);
			for (int k = 0; k < 16 * (int)1.000; k++)
			{
				ckks_encoder.encode(temps[k], levelOne, levelOneScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesOriCip[j - 1][k]);
			}
			dataDivideEncode(noises[j - 1], temps, temps + 8 * (int)1.000, batchsize * time / 8, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), empty, 8 * (int)1.000, 0);
			for (int k = 0; k < 16 * (int)1.000; k++)
			{
				ckks_encoder.encode(temps[k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesCip[j - 1][k]);
			}
		}
		//数据提供者：将16个密文的后半部分广播/前半部分发送，16×2×(partyNum - 1)个噪声密文发送，16×(partyNum - 1)个部分解密发送
		//研究者：完全解密并重新编码(将16个layer_1_2_1OriCip发送给数据提供者)
		for (int j = 0; j < 16 * (int)1.000; j++)
		{
			decryptor[0]->final_decrypt(layer_1_2_1Cip[j], pDecPlain[j], partyNum - 1, tempPlain);
			ckks_encoder.decode(tempPlain, noises[0][j]);
		}
		dataDivide(noises[0], temps, temps + 8 * (int)1.000, batchsize * time / 8, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), empty, 8 * (int)1.000, 1);
		for (int j = 0; j < 16 * (int)1.000; j++)
		{
			ckks_encoder.encode(temps[j], levelOne, levelOneScale, tempPlain);
			layer_1_2_1OriCip[j] = noisesOriCip[0][j];
			evaluator.add_plain_inplace(layer_1_2_1OriCip[j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_2_1OriCip[j], noisesOriCip[k - 1][j]);
				k++;
			}
		}
		dataDivideEncode(noises[0], temps, temps + 8 * (int)1.000, batchsize * time / 8, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), empty, 8 * (int)1.000, 1);
		for (int j = 0; j < 16 * (int)1.000; j++)
		{
			ckks_encoder.encode(temps[j], levelTwo, levelTwoScale, tempPlain);
			layer_1_2_1Cip[j] = noisesCip[0][j];
			evaluator.add_plain_inplace(layer_1_2_1Cip[j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_2_1Cip[j], noisesCip[k - 1][j]);
				k++;
			}
		}
		//研究者：正向传播8,24,40
		for (int j = 0; j < 4; j++)
		{
			evaluator.multiply_plain(wb_3Cip, whb_filterPlain[j], wb_3EncodeCip[j]);
			for (int k = pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(wb_3EncodeCip[j], -k, galois_keys, tempCip);
				evaluator.add_inplace(wb_3EncodeCip[j], tempCip);
			}
			if (j != 0)
			{
				evaluator.rotate_vector_inplace(wb_3EncodeCip[j], j * pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4, galois_keys);
			}
			evaluator.rescale_to_next_inplace(wb_3EncodeCip[j]);
		}
		parameterRotate(wb_3EncodeCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), paramCip[0], evaluator, galois_keys, one_1Plain);
		matrixMultiple_CipCip(paramCip[0], layer_1_2_1Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 8 * (int)1.000, o1_3_0Cip, evaluator, galois_keys);
		for (int j = 0; j < 8 * (int)1.000; j++)
		{
			evaluator.relinearize_inplace(o1_3_0Cip[j], relin_keys);
			evaluator.rescale_to_next_inplace(o1_3_0Cip[j]);
			evaluator.square(o1_3_0Cip[j], layer_1_3_0Cip[j]);
			evaluator.relinearize_inplace(layer_1_3_0Cip[j], relin_keys);
			evaluator.multiply_plain(o1_3_0Cip[j], yicixiang_3Plain, tempCip);
			evaluator.add_inplace(layer_1_3_0Cip[j], tempCip);
			evaluator.rescale_to_next_inplace(layer_1_3_0Cip[j]);
		}
		//研究者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < 8 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					noises[j - 1][k][m] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
				}
				ckks_encoder.encode(noises[j - 1][k], levelFour, levelFourScale, tempPlain);
				decryptor[j]->partial_decrypt(layer_1_3_0Cip[k], tempPlain, pDecPlain[k][j - 1]);
				for (int m = 0; m < slotCount; m++)
				{
					noises[j - 1][k][m] /= -pow(2, 7.75);
				}
				ckks_encoder.encode(noises[j - 1][k], levelOne, levelOneScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesOriCip[j - 1][k]);
			}
			matrixRefresh(noises[j - 1], temps, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 8 * (int)1.000, batchsize * time / 16);
			for (int k = 0; k < 8 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					temps[k][m] *= 3.735111370627234931 / pow(2, 7.75);
				}
				ckks_encoder.encode(temps[k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesCip[j - 1][k]);
			}
		}
		//研究者：将8个密文的后半部分广播/前半部分发送，8×2×(partyNum - 1)个噪声密文发送，8×(partyNum - 1)个部分解密发送
		//数据提供者：完全解密并重新编码
		for (int j = 0; j < 8 * (int)1.000; j++)
		{
			decryptor[0]->final_decrypt(layer_1_3_0Cip[j], pDecPlain[j], partyNum - 1, tempPlain);
			ckks_encoder.decode(tempPlain, noises[0][j]);
		}
		for (int j = 0; j < 8 * (int)1.000; j++)
		{
			for (int k = 0; k < slotCount; k++)
			{
				noises[0][j][k] /= pow(2, 7.75);
			}
			ckks_encoder.encode(noises[0][j], levelOne, levelOneScale, tempPlain);
			layer_1_3_0OriCip[j] = noisesOriCip[0][j];
			evaluator.add_plain_inplace(layer_1_3_0OriCip[j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_3_0OriCip[j], noisesOriCip[k - 1][j]);
				k++;
			}
		}
		matrixRefresh(noises[0], temps, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 8 * (int)1.000, batchsize * time / 16);
		for (int j = 0; j < 8 * (int)1.000; j++)
		{
			for (int k = 0; k < slotCount; k++)
			{
				temps[j][k] *= 3.735111370627234931 / pow(2, 7.75);
			}
			ckks_encoder.encode(temps[j], levelTwo, levelTwoScale, tempPlain);
			layer_1_3_0Cip[j] = noisesCip[0][j];
			evaluator.add_plain_inplace(layer_1_3_0Cip[j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_3_0Cip[j], noisesCip[k - 1][j]);
				k++;
			}
		}
		//数据提供者：计算wh_3 * layer_1_3_0
		for (int j = 0; j < 4; j++)
		{
			evaluator.multiply_plain(wh_3Cip, whb_filterPlain[j], wh_3EncodeCip[j]);
			for (int k = pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(wh_3EncodeCip[j], -k, galois_keys, tempCip);
				evaluator.add_inplace(wh_3EncodeCip[j], tempCip);
			}
			if (j != 0)
			{
				evaluator.rotate_vector_inplace(wh_3EncodeCip[j], j * pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4, galois_keys);
			}
			evaluator.rescale_to_next_inplace(wh_3EncodeCip[j]);
		}
		parameterRotate(wh_3EncodeCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), paramCip[1], evaluator, galois_keys, one_1Plain);
		matrixMultiple_CipCip(paramCip[1], layer_1_3_0Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 8 * (int)1.000, layer_1_3_1Cip, evaluator, galois_keys);
		for (int j = 0; j < 8 * (int)1.000; j++)
		{
			evaluator.relinearize_inplace(layer_1_3_1Cip[j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_3_1Cip[j]);
		}
		//数据提供者：将8个密文发送
		//研究者：正向传播16,32,48
		matrixMultiple_CipCip(paramCip[0], layer_1_2_1Cip + 8 * (int)1.000, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 8 * (int)1.000, o1_3_1Cip, evaluator, galois_keys);
		for (int j = 0; j < 8 * (int)1.000; j++)
		{
			evaluator.relinearize_inplace(o1_3_1Cip[j], relin_keys);
			evaluator.rescale_to_next_inplace(o1_3_1Cip[j]);
			evaluator.add_inplace(o1_3_1Cip[j], layer_1_3_1Cip[j]);
			evaluator.square(o1_3_1Cip[j], layer_1_3_1Cip[j]);
			evaluator.relinearize_inplace(layer_1_3_1Cip[j], relin_keys);
			evaluator.multiply_plain(o1_3_1Cip[j], yicixiang_3Plain, tempCip);
			evaluator.add_inplace(layer_1_3_1Cip[j], tempCip);
			evaluator.rescale_to_next_inplace(layer_1_3_1Cip[j]);
		}
		//研究者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < 8 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					noises[j - 1][k][m] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
				}
				ckks_encoder.encode(noises[j - 1][k], levelFour, levelFourScale, tempPlain);
				decryptor[j]->partial_decrypt(layer_1_3_1Cip[k], tempPlain, pDecPlain[k][j - 1]);
				for (int m = 0; m < slotCount; m++)
				{
					noises[j - 1][k][m] *= -1;
				}
			}
			dataDivideThree(noises[j - 1], temps, temps + 4 * (int)1.000, temps + 8 * (int)1.000, batchsize * time / 16, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), empty, 4 * (int)1.000, 0);
			for (int k = 0; k < 12 * (int)1.000; k++)
			{
				ckks_encoder.encode(temps[k], levelOne, levelOneScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesOriCip[j - 1][k]);
			}
			dataDivideEncodeThree(noises[j - 1], temps, temps + 4 * (int)1.000, temps + 8 * (int)1.000, batchsize * time / 16, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), empty, 4 * (int)1.000, 0);
			for (int k = 0; k < 12 * (int)1.000; k++)
			{
				ckks_encoder.encode(temps[k], levelOne, levelOneScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesCip[j - 1][k]);
			}
		}
		//研究者：将8个密文的后半部分广播/前半部分发送，12×2×(partyNum - 1)个噪声密文发送，8×(partyNum - 1)个部分解密发送
		//数据提供者：完全解密并重新编码
		for (int j = 0; j < 8 * (int)1.000; j++)
		{
			decryptor[0]->final_decrypt(layer_1_3_1Cip[j], pDecPlain[j], partyNum - 1, tempPlain);
			ckks_encoder.decode(tempPlain, noises[0][j]);
		}
		dataDivideThree(noises[0], temps, temps + 4 * (int)1.000, temps + 8 * (int)1.000, batchsize * time / 16, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), empty, 4 * (int)1.000, 1);
		for (int j = 0; j < 12 * (int)1.000; j++)
		{
			ckks_encoder.encode(temps[j], levelOne, levelOneScale, tempPlain);
			layer_1_3_1OriCip[j] = noisesOriCip[0][j];
			evaluator.add_plain_inplace(layer_1_3_1OriCip[j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_3_1OriCip[j], noisesOriCip[k - 1][j]);
				k++;
			}
		}
		dataDivideEncodeThree(noises[0], temps, temps + 4 * (int)1.000, temps + 8 * (int)1.000, batchsize * time / 16, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), empty, 4 * (int)1.000, 1);
		for (int j = 0; j < 12 * (int)1.000; j++)
		{
			ckks_encoder.encode(temps[j], levelOne, levelOneScale, tempPlain);
			layer_1_3_1Cip[j] = noisesCip[0][j];
			evaluator.add_plain_inplace(layer_1_3_1Cip[j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_3_1Cip[j], noisesCip[k - 1][j]);
				k++;
			}
		}
		//数据提供者：正向传播16
		for (int j = 0; j < 4; j++)
		{
			evaluator.multiply_plain(wb_4Cip, whb_filterPlain[j], wb_4EncodeCip[j]);
			for (int k = pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(wb_4EncodeCip[j], -k, galois_keys, tempCip);
				evaluator.add_inplace(wb_4EncodeCip[j], tempCip);
			}
			if (j != 0)
			{
				evaluator.rotate_vector_inplace(wb_4EncodeCip[j], j * pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4, galois_keys);
			}
			evaluator.rescale_to_next_inplace(wb_4EncodeCip[j]);
		}
		parameterRotate(wb_4EncodeCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), paramCip[0], evaluator, galois_keys, yicixiang_3Plain);
		matrixMultiple_CipCip(paramCip[0], layer_1_3_1Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 4 * (int)1.000, o1_4_0Cip, evaluator, galois_keys);
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			evaluator.relinearize_inplace(o1_4_0Cip[j], relin_keys);
			evaluator.rescale_to_next_inplace(o1_4_0Cip[j]);
			evaluator.square(o1_4_0Cip[j], layer_1_4_0Cip[j]);
			evaluator.relinearize_inplace(layer_1_4_0Cip[j], relin_keys);
			evaluator.multiply_plain(o1_4_0Cip[j], yicixiang_2Plain, tempCip);
			evaluator.add_inplace(layer_1_4_0Cip[j], tempCip);
			evaluator.rescale_to_next_inplace(layer_1_4_0Cip[j]);
		}
		//数据提供者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < 4 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					noises[j - 1][k][m] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
				}
				ckks_encoder.encode(noises[j - 1][k], levelThree, levelThreeScale, tempPlain);
				decryptor[j]->partial_decrypt(layer_1_4_0Cip[k], tempPlain, pDecPlain[k][j - 1]);
				for (int m = 0; m < slotCount; m++)
				{
					noises[j - 1][k][m] /= -pow(2, 7.75);
				}
				ckks_encoder.encode(noises[j - 1][k], levelOne, levelOneScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesOriCip[j - 1][k]);
			}
			matrixRefresh(noises[j - 1], temps, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 4 * (int)1.000, batchsize);
			for (int k = 0; k < 4 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					temps[k][m] *= 3.735111370627234931 / pow(2, 7.75);
				}
				ckks_encoder.encode(temps[k], levelOne, levelOneScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesCip[j - 1][k]);
			}
		}
		//数据提供者：将4个密文的后半部分广播/前半部分发送，4×2×(partyNum - 1)个噪声密文发送，4×(partyNum - 1)个部分解密发送
		//研究者：完全解密并重新编码(将4个layer_1_4_0OriCip发送给数据提供者)
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			decryptor[0]->final_decrypt(layer_1_4_0Cip[j], pDecPlain[j], partyNum - 1, tempPlain);
			ckks_encoder.decode(tempPlain, noises[0][j]);
		}
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			for (int k = 0; k < slotCount; k++)
			{
				noises[0][j][k] /= pow(2, 7.75);
			}
			ckks_encoder.encode(noises[0][j], levelOne, levelOneScale, tempPlain);
			layer_1_4_0OriCip[j] = noisesOriCip[0][j];
			evaluator.add_plain_inplace(layer_1_4_0OriCip[j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_4_0OriCip[j], noisesOriCip[k - 1][j]);
				k++;
			}
		}
		matrixRefresh(noises[0], temps, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 4 * (int)1.000, batchsize);
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			for (int k = 0; k < slotCount; k++)
			{
				temps[j][k] *= 3.735111370627234931 / pow(2, 7.75);
			}
			ckks_encoder.encode(temps[j], levelOne, levelOneScale, tempPlain);
			layer_1_4_0Cip[j] = noisesCip[0][j];
			evaluator.add_plain_inplace(layer_1_4_0Cip[j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_4_0Cip[j], noisesCip[k - 1][j]);
				k++;
			}
		}
		//研究者：计算wh_4 * layer_1_4_0
		for (int j = 0; j < 4; j++)
		{
			evaluator.multiply_plain(wh_4Cip, whb_filterPlain[j], wh_4EncodeCip[j]);
			for (int k = pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(wh_4EncodeCip[j], -k, galois_keys, tempCip);
				evaluator.add_inplace(wh_4EncodeCip[j], tempCip);
			}
			if (j != 0)
			{
				evaluator.rotate_vector_inplace(wh_4EncodeCip[j], j * pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4, galois_keys);
			}
			evaluator.rescale_to_next_inplace(wh_4EncodeCip[j]);
		}
		parameterRotate(wh_4EncodeCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), paramCip[1], evaluator, galois_keys, yicixiang_3Plain);
		matrixMultiple_CipCip(paramCip[1], layer_1_4_0Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 4 * (int)1.000, layer_1_4_1Cip, evaluator, galois_keys);
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			evaluator.relinearize_inplace(layer_1_4_1Cip[j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_4_1Cip[j]);
		}
		//研究者：将4个密文发送
		//数据提供者：正向传播32
		matrixMultiple_CipCip(paramCip[0], layer_1_3_1Cip + 4 * (int)1.000, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 4 * (int)1.000, o1_4_1Cip, evaluator, galois_keys);
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			evaluator.relinearize_inplace(o1_4_1Cip[j], relin_keys);
			evaluator.rescale_to_next_inplace(o1_4_1Cip[j]);
			evaluator.add_inplace(o1_4_1Cip[j], layer_1_4_1Cip[j]);
			evaluator.square(o1_4_1Cip[j], layer_1_4_1Cip[j]);
			evaluator.relinearize_inplace(layer_1_4_1Cip[j], relin_keys);
			evaluator.multiply_plain(o1_4_1Cip[j], yicixiang_2Plain, tempCip);
			evaluator.add_inplace(layer_1_4_1Cip[j], tempCip);
			evaluator.rescale_to_next_inplace(layer_1_4_1Cip[j]);
		}
		//数据提供者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < 4 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					noises[j - 1][k][m] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
				}
				ckks_encoder.encode(noises[j - 1][k], levelThree, levelThreeScale, tempPlain);
				decryptor[j]->partial_decrypt(layer_1_4_1Cip[k], tempPlain, pDecPlain[k][j - 1]);
				for (int m = 0; m < slotCount; m++)
				{
					noises[j - 1][k][m] /= -pow(2, 7.75);
				}
				ckks_encoder.encode(noises[j - 1][k], levelOne, levelOneScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesOriCip[j - 1][k]);
			}
			matrixRefresh(noises[j - 1], temps, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 4 * (int)1.000, batchsize);
			for (int k = 0; k < 4 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					temps[k][m] *= 3.735111370627234931 / pow(2, 7.75);
				}
				ckks_encoder.encode(temps[k], levelOne, levelOneScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesCip[j - 1][k]);
			}
		}
		//数据提供者：将4个密文的后半部分广播/前半部分发送，4×2×(partyNum - 1)个噪声密文发送，4×(partyNum - 1)个部分解密发送
		//研究者：完全解密并重新编码(将4个layer_1_4_1OriCip发送给数据提供者)
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			decryptor[0]->final_decrypt(layer_1_4_1Cip[j], pDecPlain[j], partyNum - 1, tempPlain);
			ckks_encoder.decode(tempPlain, noises[0][j]);
		}
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			for (int k = 0; k < slotCount; k++)
			{
				noises[0][j][k] /= pow(2, 7.75);
			}
			ckks_encoder.encode(noises[0][j], levelOne, levelOneScale, tempPlain);
			layer_1_4_1OriCip[j] = noisesOriCip[0][j];
			evaluator.add_plain_inplace(layer_1_4_1OriCip[j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_4_1OriCip[j], noisesOriCip[k - 1][j]);
				k++;
			}
		}
		matrixRefresh(noises[0], temps, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 4 * (int)1.000, batchsize);
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			for (int k = 0; k < slotCount; k++)
			{
				temps[j][k] *= 3.735111370627234931 / pow(2, 7.75);
			}
			ckks_encoder.encode(temps[j], levelOne, levelOneScale, tempPlain);
			layer_1_4_1Cip[j] = noisesCip[0][j];
			evaluator.add_plain_inplace(layer_1_4_1Cip[j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_4_1Cip[j], noisesCip[k - 1][j]);
				k++;
			}
		}
		//研究者：计算wh_4 * layer_1_4_1
		matrixMultiple_CipCip(paramCip[1], layer_1_4_1Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 4 * (int)1.000, layer_1_4_2Cip, evaluator, galois_keys);
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			evaluator.relinearize_inplace(layer_1_4_2Cip[j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_4_2Cip[j]);
		}
		//研究者：将4个密文发送
		//数据提供者：正向传播48
		matrixMultiple_CipCip(paramCip[0], layer_1_3_1Cip + 8 * (int)1.000, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 4 * (int)1.000, o1_4_2Cip, evaluator, galois_keys);
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			evaluator.relinearize_inplace(o1_4_2Cip[j], relin_keys);
			evaluator.rescale_to_next_inplace(o1_4_2Cip[j]);
			evaluator.add_inplace(o1_4_2Cip[j], layer_1_4_2Cip[j]);
			evaluator.square(o1_4_2Cip[j], layer_1_4_2Cip[j]);
			evaluator.relinearize_inplace(layer_1_4_2Cip[j], relin_keys);
			evaluator.multiply_plain(o1_4_2Cip[j], yicixiang_2Plain, tempCip);
			evaluator.add_inplace(layer_1_4_2Cip[j], tempCip);
			evaluator.rescale_to_next_inplace(layer_1_4_2Cip[j]);
		}
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 1 * (int)1.000; k++)
			{
				evaluator.add_plain_inplace(layer_1_4_2Cip[2 + j + k * 4], one_3Plain);
			}
		}
		//数据提供者：计算w1 * layer_1_4_2
		for (int j = 0; j < 2; j++)
		{
			evaluator.multiply_plain(wb1_Cip, wb1_filterPlain[j], wb1_EncodeCip[j]);
			evaluator.rotate_vector(wb1_EncodeCip[j], (-slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + pow(2, ceil(log2(hidenode * hidenode / 4) / 2))) / 2, galois_keys, tempCip);
			evaluator.add_inplace(wb1_EncodeCip[j], tempCip);
			if (j != 0)
			{
				evaluator.rotate_vector_inplace(wb1_EncodeCip[j], (long long)slotCount * j / pow(2, 1 + ceil(log2(hidenode * hidenode / 4) / 2)), galois_keys);
			}
			evaluator.rescale_to_next_inplace(wb1_EncodeCip[j]);
			evaluator.multiply_plain_inplace(wb1_EncodeCip[j], one_1Plain);
			evaluator.rescale_to_next_inplace(wb1_EncodeCip[j]);
			evaluator.multiply_plain(wb1_EncodeCip[j], one_2Plain, paramCip[0][j]);
			evaluator.rescale_to_next_inplace(paramCip[0][j]);
		}
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			evaluator.multiply(layer_1_4_2Cip[j], paramCip[0][j % 4 / 2], paramCip[1][j]);
		}
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 1 * (int)1.000; k++)
			{
				evaluator.add(paramCip[1][j + k * 4], paramCip[1][j + 2 + k * 4], layer_2Cip[j + k * 2]);
				evaluator.relinearize_inplace(layer_2Cip[j + k * 2], relin_keys);
			}
		}
		for (int j = 0; j < 2 * (int)1.000; j++)
		{
			for (int k = slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k < slotCount; k *= 2)
			{
				evaluator.rotate_vector(layer_2Cip[j], k, galois_keys, tempCip);
				evaluator.add_inplace(layer_2Cip[j], tempCip);
			}
			evaluator.rescale_to_next_inplace(layer_2Cip[j]);
			for (int k = 0; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k++)
			{
				layer_2_error[j][k] = (((rand() << 15) + rand()) % 11000 / 1000.0 - 5.5) * pow(2, 15.5);
				for (int m = 1; m < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); m++)
				{
					layer_2_error[j][k + (long long)m * slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2))] = layer_2_error[j][k];
				}
			}
			ckks_encoder.encode(layer_2_error[j], levelFour, levelFourScale, layer_2_errorPlain[j]);
			evaluator.add_plain_inplace(layer_2Cip[j], layer_2_errorPlain[j]);
			for (int k = 0; k < slotCount / pow(2, 1 + ceil(log2(hidenode * hidenode / 4) / 2)); k++)
			{
				for (int m = 0; m < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); m++)
				{
					layer_2_error[j][k + (long long)m * slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + slotCount / pow(2, 1 + ceil(log2(hidenode * hidenode / 4) / 2))] = layer_2_error[j][k];
				}
			}
			for (int k = 0; k < slotCount; k++)
			{
				layer_2_error[j][k] /= pow(2, 31.0 / 3.0);
			}
			ckks_encoder.encode(layer_2_error[j], levelZeroScale, layer_2_errorPlain[j]);
			xieLv[j] = empty;
			jieJu[j] = empty;
			for (int k = 0; k < slotCount / pow(2, 1 + ceil(log2(hidenode * hidenode / 4) / 2)); k++)
			{
				for (int m = 0; m < pow(2, 1 + ceil(log2(hidenode * hidenode / 4) / 2)); m++)
				{
					xieLv[j][k + (long long)m * slotCount / (long long)pow(2, 1 + ceil(log2(hidenode * hidenode / 4) / 2))] = xieLvKu[min(max(m - (int)round(layer_2_error[j][k] / pow(2, 15.5 / 3.0)), 0), 31)];
					jieJu[j][k + (long long)m * slotCount / (long long)pow(2, 1 + ceil(log2(hidenode * hidenode / 4) / 2))] = jieJuKu[min(max(m - (int)round(layer_2_error[j][k] / pow(2, 15.5 / 3.0)), 0), 31)] * pow(2, 15.5 / 3.0);
				}
				ckks_encoder.encode(xieLv[j], levelZeroScale, xieLvPlain[j]);
				ckks_encoder.encode(jieJu[j], levelOne, levelOneScale, jieJuPlain[j]);
			}
		}
		//数据提供者：增加噪声并部分解密
		temp = empty;
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < 2 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					noises[j - 1][k][m] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
					temp[(long long)k * slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m % (long long)(slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))] += noises[j - 1][k][m];
				}
				ckks_encoder.encode(noises[j - 1][k], levelFour, levelFourScale, tempPlain);
				decryptor[j]->partial_decrypt(layer_2Cip[k], tempPlain, pDecPlain[k][j - 1]);
			}
		}
		//数据提供者：将2个密文的后半部分广播/前半部分发送，2×(partyNum - 1)个部分解密发送
		//研究者：完全解密并重新编码
		for (int j = 0; j < 2 * (int)1.000; j++)
		{
			decryptor[0]->final_decrypt(layer_2Cip[j], pDecPlain[j], partyNum - 1, tempPlain);
			ckks_encoder.decode(tempPlain, noises[0][j]);
			temps[j] = empty;
		}
		noises[1][0] = empty;
		for (int j = 0; j < batchsize; j++)
		{
			for (int k = 0; k < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k++)
			{
				noises[1][0][j] += noises[0][j / (long long)(slotCount / pow(2, 1 + ceil(log2(hidenode * hidenode / 4) / 2)))][j % (long long)(slotCount / pow(2, 1 + ceil(log2(hidenode * hidenode / 4) / 2))) + (long long)k * slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2))];
			}
			noises[1][0][j] -= temp[j + j / (long long)(slotCount / pow(2, 1 + ceil(log2(hidenode * hidenode / 4) / 2))) * slotCount / pow(2, 1 + ceil(log2(hidenode * hidenode / 4) / 2))];
			noises[1][0][j] /= pow(2, ceil(log2(hidenode * hidenode / 4) / 2));
			for (int k = 0; k < pow(2, 1 + ceil(log2(hidenode * hidenode / 4) / 2)); k++)
			{
				temps[j / (long long)(slotCount / pow(2, 1 + ceil(log2(hidenode * hidenode / 4) / 2)))][j % (long long)(slotCount / pow(2, 1 + ceil(log2(hidenode * hidenode / 4) / 2))) + (long long)k * slotCount / pow(2, 1 + ceil(log2(hidenode * hidenode / 4) / 2))] = noises[1][0][j];
			}
		}
		for (int j = 0; j < 2 * (int)1.000; j++)
		{
			for (int k = 0; k < slotCount; k++)
			{
				temps[j][k] /= pow(2, 31.0 / 3.0);
			}
			ckks_encoder.encode(temps[j], levelZeroScale, tempPlain);
			encryptor.encrypt(tempPlain, layer_2Cip[j]);
			temps[j] = empty;
		}
		for (int j = 0; j < batchsize; j++)
		{
			temps[j / (long long)(slotCount / pow(2, 1 + ceil(log2(hidenode * hidenode / 4) / 2)))][j % (long long)(slotCount / pow(2, 1 + ceil(log2(hidenode * hidenode / 4) / 2))) + min(max(0.0, floor(noises[1][0][j] / pow(2, 15.5)) + 16), 31.0) * slotCount / pow(2, 1 + ceil(log2(hidenode * hidenode / 4) / 2))] = 1;
		}
		for (int j = 0; j < 2 * (int)1.000; j++)
		{
			ckks_encoder.encode(temps[j], levelOne, levelOneScale, tempPlain);
			encryptor.encrypt(tempPlain, noisesCip[0][j]);
		}
		//研究者：将2个输出密文，2个位置指示密文发送
		//数据提供者：计算layer_2和layer_2_delta，反向传播48_4
		for (int j = 0; j < 2 * (int)1.000; j++)
		{
			evaluator.sub_plain_inplace(layer_2Cip[j], layer_2_errorPlain[j]);
			evaluator.multiply_plain_inplace(layer_2Cip[j], xieLvPlain[j]);
			evaluator.rescale_to_next_inplace(layer_2Cip[j]);
			evaluator.add_plain_inplace(layer_2Cip[j], jieJuPlain[j]);
			evaluator.multiply_inplace(layer_2Cip[j], noisesCip[0][j]);
			evaluator.relinearize_inplace(layer_2Cip[j], relin_keys);
			for (int k = slotCount / pow(2, 1 + ceil(log2(hidenode * hidenode / 4) / 2)); k < slotCount; k *= 2)
			{
				evaluator.rotate_vector(layer_2Cip[j], k, galois_keys, tempCip);
				evaluator.add_inplace(layer_2Cip[j], tempCip);
			}
			evaluator.rescale_to_next_inplace(layer_2Cip[j]);
		}
		for (int j = 0; j < 2 * (int)1.000; j++)
		{
			y[j] = empty;
		}
		Pos = 0;
		for (int j = 0; j < batchsize; j++)
		{
			for (int k = 0; k < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k++)
			{
				y[Pos + j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) % 2][j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2 + 1)) * pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j % (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + (long long)k * slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2))] = outcome[(i * batchsize + j) % row] * pow(2, 15.5 / 3.0);
			}
			if (j % (int)((slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) - pow(2, ceil(log2(hidenode * hidenode / 4) / 2))) * 2) == (slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) - pow(2, ceil(log2(hidenode * hidenode / 4) / 2))) * 2 - 1)
			{
				Pos += 2;
			}
		}
		for (int j = 0; j < 2 * (int)1.000; j++)
		{
			ckks_encoder.encode(y[j], levelTwo, levelTwoScale, tempPlain);
			evaluator.sub_plain(layer_2Cip[j], tempPlain, layer_2_deltaCip[j]);
			evaluator.negate_inplace(layer_2_deltaCip[j]);
		}
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			evaluator.multiply(layer_2_deltaCip[j % 2 + j / 4 * 2], wb1_EncodeCip[j % 4 / 2], layer_1_4_deltaCip[2][j]);
			evaluator.relinearize_inplace(layer_1_4_deltaCip[2][j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_4_deltaCip[2][j]);
			evaluator.multiply_plain_inplace(o1_4_2Cip[j], one_2Plain);
			evaluator.rescale_to_next_inplace(o1_4_2Cip[j]);
			evaluator.add_plain_inplace(o1_4_2Cip[j], half_3Plain);
			evaluator.multiply_inplace(layer_1_4_deltaCip[2][j], o1_4_2Cip[j]);
			evaluator.relinearize_inplace(layer_1_4_deltaCip[2][j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_4_deltaCip[2][j]);
		}
		//数据提供者：将4个密文的后半部分广播
		//研究者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < 4 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[2][j - 1][k][m] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
				}
				ckks_encoder.encode(backnoises[2][j - 1][k], levelFour, levelFourScale, tempPlain);
				decryptor[j]->partial_decrypt(layer_1_4_deltaCip[2][k], tempPlain, pDecPlain[k][j - 1]);
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[2][j - 1][k][m] *= -27.902113901977722956 / pow(2, 77.5 / 6.0);
				}
				ckks_encoder.encode(backnoises[2][j - 1][k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesOriCip[j - 1][k]);
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[2][j - 1][k][m] /= 5.282245157314995551;
				}
			}
			matrixRefresh(backnoises[2][j - 1], temps, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 4 * (int)1.000, batchsize);
			for (int k = 0; k < 4 * (int)1.000; k++)
			{
				ckks_encoder.encode(temps[k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesCip[j - 1][k]);
			}
		}
		//研究者：4×2×(partyNum - 1)个噪声密文发送，4×(partyNum - 1)个部分解密发送
		//数据提供者：完全解密并重新编码
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			decryptor[0]->final_decrypt(layer_1_4_deltaCip[2][j], pDecPlain[j], partyNum - 1, tempPlain);
			ckks_encoder.decode(tempPlain, noises[2][j]);
		}
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			for (int k = 0; k < slotCount; k++)
			{
				noises[2][j][k] *= 27.902113901977722956 / pow(2, 77.5 / 6.0);
			}
			ckks_encoder.encode(noises[2][j], levelTwo, levelTwoScale, tempPlain);
			layer_1_4_deltaOriCip[2][j][0] = noisesOriCip[0][j];
			evaluator.add_plain_inplace(layer_1_4_deltaOriCip[2][j][0], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_4_deltaOriCip[2][j][0], noisesOriCip[k - 1][j]);
				k++;
			}
			for (int k = 0; k < slotCount; k++)
			{
				noises[2][j][k] /= 5.282245157314995551;
			}
		}
		matrixRefresh(noises[2], temps, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 4 * (int)1.000, batchsize);
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			ckks_encoder.encode(temps[j], levelTwo, levelTwoScale, tempPlain);
			layer_1_4_deltaCip[2][j] = noisesCip[0][j];
			evaluator.add_plain_inplace(layer_1_4_deltaCip[2][j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_4_deltaCip[2][j], noisesCip[k - 1][j]);
				k++;
			}
		}
		//数据提供者：反向传播32_4
		for (int j = 0; j < 4; j++)
		{
			evaluator.multiply_plain(wh_4TranCip, whb_filterPlain[j], wh_4EncodeCip[j]);
			for (int k = pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(wh_4EncodeCip[j], -k, galois_keys, tempCip);
				evaluator.add_inplace(wh_4EncodeCip[j], tempCip);
			}
			if (j != 0)
			{
				evaluator.rotate_vector_inplace(wh_4EncodeCip[j], j * pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4, galois_keys);
			}
			evaluator.rescale_to_next_inplace(wh_4EncodeCip[j]);
		}
		parameterRotate(wh_4EncodeCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), paramCip[1], evaluator, galois_keys, one_1Plain);
		matrixMultiple_CipCip(paramCip[1], layer_1_4_deltaCip[2], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 4 * (int)1.000, layer_1_4_deltaCip[1], evaluator, galois_keys);
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			evaluator.relinearize_inplace(layer_1_4_deltaCip[1][j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_4_deltaCip[1][j]);
			evaluator.multiply_plain_inplace(o1_4_1Cip[j], one_2Plain);
			evaluator.rescale_to_next_inplace(o1_4_1Cip[j]);
			evaluator.add_plain_inplace(o1_4_1Cip[j], half_3Plain);
			evaluator.multiply_inplace(layer_1_4_deltaCip[1][j], o1_4_1Cip[j]);
			evaluator.relinearize_inplace(layer_1_4_deltaCip[1][j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_4_deltaCip[1][j]);
		}
		//数据提供者：将4个密文的后半部分广播
		//研究者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < 4 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[1][j - 1][k][m] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
				}
				ckks_encoder.encode(backnoises[1][j - 1][k], levelFour, levelFourScale, tempPlain);
				decryptor[j]->partial_decrypt(layer_1_4_deltaCip[1][k], tempPlain, pDecPlain[k][j - 1]);
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[1][j - 1][k][m] *= -27.902113901977722956 / pow(2, 15.5);
				}
				ckks_encoder.encode(backnoises[1][j - 1][k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesOriCip[j - 1][k]);
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[1][j - 1][k][m] /= 5.282245157314995551;
				}
			}
			matrixRefresh(backnoises[1][j - 1], temps, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 4 * (int)1.000, batchsize);
			for (int k = 0; k < 4 * (int)1.000; k++)
			{
				ckks_encoder.encode(temps[k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesCip[j - 1][k]);
			}
		}
		//研究者：将4×2×(partyNum - 1)个噪声密文发送，4×(partyNum - 1)个部分解密发送
		//数据提供者：完全解密并重新编码
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			decryptor[0]->final_decrypt(layer_1_4_deltaCip[1][j], pDecPlain[j], partyNum - 1, tempPlain);
			ckks_encoder.decode(tempPlain, noises[1][j]);
		}
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			for (int k = 0; k < slotCount; k++)
			{
				noises[1][j][k] *= 27.902113901977722956 / pow(2, 15.5);
			}
			ckks_encoder.encode(noises[1][j], levelTwo, levelTwoScale, tempPlain);
			layer_1_4_deltaOriCip[1][j][0] = noisesOriCip[0][j];
			evaluator.add_plain_inplace(layer_1_4_deltaOriCip[1][j][0], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_4_deltaOriCip[1][j][0], noisesOriCip[k - 1][j]);
				k++;
			}
			for (int k = 0; k < slotCount; k++)
			{
				noises[1][j][k] /= 5.282245157314995551;
			}
		}
		matrixRefresh(noises[1], temps, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 4 * (int)1.000, batchsize);
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			ckks_encoder.encode(temps[j], levelTwo, levelTwoScale, tempPlain);
			layer_1_4_deltaCip[1][j] = noisesCip[0][j];
			evaluator.add_plain_inplace(layer_1_4_deltaCip[1][j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_4_deltaCip[1][j], noisesCip[k - 1][j]);
				k++;
			}
		}
		//数据提供者：反向传播16_4
		matrixMultiple_CipCip(paramCip[1], layer_1_4_deltaCip[1], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 4 * (int)1.000, layer_1_4_deltaCip[0], evaluator, galois_keys);
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			evaluator.relinearize_inplace(layer_1_4_deltaCip[0][j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_4_deltaCip[0][j]);
			evaluator.multiply_plain_inplace(o1_4_0Cip[j], one_2Plain);
			evaluator.rescale_to_next_inplace(o1_4_0Cip[j]);
			evaluator.add_plain_inplace(o1_4_0Cip[j], half_3Plain);
			evaluator.multiply_inplace(layer_1_4_deltaCip[0][j], o1_4_0Cip[j]);
			evaluator.relinearize_inplace(layer_1_4_deltaCip[0][j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_4_deltaCip[0][j]);
		}
		//数据提供者：将4个密文的后半部分广播
		//研究者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < 4 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[0][j - 1][k][m] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
				}
				ckks_encoder.encode(backnoises[0][j - 1][k], levelFour, levelFourScale, tempPlain);
				decryptor[j]->partial_decrypt(layer_1_4_deltaCip[0][k], tempPlain, pDecPlain[k][j - 1]);
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[0][j - 1][k][m] *= -27.902113901977722956 / pow(2, 15.5);
				}
				ckks_encoder.encode(backnoises[0][j - 1][k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesOriCip[j - 1][k]);
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[0][j - 1][k][m] /= 5.282245157314995551;
				}
			}
			dataCombineEncodeThree(temps, backnoises[0][j - 1], backnoises[1][j - 1], backnoises[2][j - 1], batchsize * time / 16, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), empty, 8 * (int)1.000);
			for (int k = 0; k < 8 * (int)1.000; k++)
			{
				ckks_encoder.encode(temps[k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesCip[j - 1][k]);
			}
		}
		//研究者：将(8 + 4)×(partyNum - 1)个噪声密文发送，4×(partyNum - 1)个部分解密发送
		//数据提供者：完全解密并重新编码
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			decryptor[0]->final_decrypt(layer_1_4_deltaCip[0][j], pDecPlain[j], partyNum - 1, tempPlain);
			ckks_encoder.decode(tempPlain, noises[0][j]);
		}
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			for (int k = 0; k < slotCount; k++)
			{
				noises[0][j][k] *= 27.902113901977722956 / pow(2, 15.5);
			}
			ckks_encoder.encode(noises[0][j], levelTwo, levelTwoScale, tempPlain);
			layer_1_4_deltaOriCip[0][j][0] = noisesOriCip[0][j];
			evaluator.add_plain_inplace(layer_1_4_deltaOriCip[0][j][0], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_4_deltaOriCip[0][j][0], noisesOriCip[k - 1][j]);
				k++;
			}
			for (int k = 0; k < slotCount; k++)
			{
				noises[0][j][k] /= 5.282245157314995551;
			}
		}
		dataCombineEncodeThree(temps, noises[0], noises[1], noises[2], batchsize * time / 16, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), empty, 8 * (int)1.000);
		for (int j = 0; j < 8 * (int)1.000; j++)
		{
			ckks_encoder.encode(temps[j], levelTwo, levelTwoScale, tempPlain);
			evaluator.add_plain_inplace(noisesCip[0][j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(noisesCip[0][j], noisesCip[k - 1][j]);
				k++;
			}
		}
		//数据提供者：将8个密文发送
		//研究者：反向传播16倍数_3
		for (int j = 0; j < 4; j++)
		{
			evaluator.multiply_plain(wb_4TranCip, whb_filterPlain[j], wb_4EncodeCip[j]);
			for (int k = pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(wb_4EncodeCip[j], -k, galois_keys, tempCip);
				evaluator.add_inplace(wb_4EncodeCip[j], tempCip);
			}
			if (j != 0)
			{
				evaluator.rotate_vector_inplace(wb_4EncodeCip[j], j * pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4, galois_keys);
			}
			evaluator.rescale_to_next_inplace(wb_4EncodeCip[j]);
		}
		parameterRotate(wb_4EncodeCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), paramCip[0], evaluator, galois_keys, one_1Plain);
		matrixMultiple_CipCip(paramCip[0], noisesCip[0], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 8 * (int)1.000, layer_1_3_deltaCip[1], evaluator, galois_keys);
		for (int j = 0; j < 8 * (int)1.000; j++)
		{
			evaluator.relinearize_inplace(layer_1_3_deltaCip[1][j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_3_deltaCip[1][j]);
			evaluator.add_plain_inplace(o1_3_1Cip[j], half_3Plain);
			evaluator.multiply_inplace(layer_1_3_deltaCip[1][j], o1_3_1Cip[j]);
			evaluator.relinearize_inplace(layer_1_3_deltaCip[1][j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_3_deltaCip[1][j]);
		}
		//研究者：将8个密文的后半部分广播
		//数据提供者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < 8 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[1][j - 1][k][m] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
				}
				ckks_encoder.encode(backnoises[1][j - 1][k], levelFour, levelFourScale, tempPlain);
				decryptor[j]->partial_decrypt(layer_1_3_deltaCip[1][k], tempPlain, pDecPlain[k][j - 1]);
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[1][j - 1][k][m] *= -27.902113901977722956 / pow(2, 15.5);
				}
				ckks_encoder.encode(backnoises[1][j - 1][k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesOriCip[j - 1][k]);
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[1][j - 1][k][m] /= 5.282245157314995551;
				}
			}
			matrixRefresh(backnoises[1][j - 1], temps, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 8 * (int)1.000, batchsize * time / 16);
			for (int k = 0; k < 8 * (int)1.000; k++)
			{
				ckks_encoder.encode(temps[k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesCip[j - 1][k]);
			}
		}
		//数据提供者：将8×2×(partyNum - 1)个噪声密文发送，8×(partyNum - 1)个部分解密发送
		//研究者：完全解密并重新编码(将8个layer_1_3_deltaOriCip[1]发送给数据提供者)
		for (int j = 0; j < 8 * (int)1.000; j++)
		{
			decryptor[0]->final_decrypt(layer_1_3_deltaCip[1][j], pDecPlain[j], partyNum - 1, tempPlain);
			ckks_encoder.decode(tempPlain, noises[1][j]);
		}
		for (int j = 0; j < 8 * (int)1.000; j++)
		{
			for (int k = 0; k < slotCount; k++)
			{
				noises[1][j][k] *= 27.902113901977722956 / pow(2, 15.5);
			}
			ckks_encoder.encode(noises[1][j], levelTwo, levelTwoScale, tempPlain);
			layer_1_3_deltaOriCip[1][j][0] = noisesOriCip[0][j];
			evaluator.add_plain_inplace(layer_1_3_deltaOriCip[1][j][0], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_3_deltaOriCip[1][j][0], noisesOriCip[k - 1][j]);
				k++;
			}
			for (int k = 0; k < slotCount; k++)
			{
				noises[1][j][k] /= 5.282245157314995551;
			}
		}
		matrixRefresh(noises[1], temps, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 8 * (int)1.000, batchsize * time / 16);
		for (int j = 0; j < 8 * (int)1.000; j++)
		{
			ckks_encoder.encode(temps[j], levelTwo, levelTwoScale, tempPlain);
			layer_1_3_deltaCip[1][j] = noisesCip[0][j];
			evaluator.add_plain_inplace(layer_1_3_deltaCip[1][j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_3_deltaCip[1][j], noisesCip[k - 1][j]);
				k++;
			}
		}
		//研究者：反向传播16余8_3
		for (int j = 0; j < 4; j++)
		{
			evaluator.multiply_plain(wh_3TranCip, whb_filterPlain[j], wh_3EncodeCip[j]);
			for (int k = pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(wh_3EncodeCip[j], -k, galois_keys, tempCip);
				evaluator.add_inplace(wh_3EncodeCip[j], tempCip);
			}
			if (j != 0)
			{
				evaluator.rotate_vector_inplace(wh_3EncodeCip[j], j * pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4, galois_keys);
			}
			evaluator.rescale_to_next_inplace(wh_3EncodeCip[j]);
		}
		parameterRotate(wh_3EncodeCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), paramCip[1], evaluator, galois_keys, one_1Plain);
		matrixMultiple_CipCip(paramCip[1], layer_1_3_deltaCip[1], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 8 * (int)1.000, layer_1_3_deltaCip[0], evaluator, galois_keys);
		for (int j = 0; j < 8 * (int)1.000; j++)
		{
			evaluator.relinearize_inplace(layer_1_3_deltaCip[0][j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_3_deltaCip[0][j]);
			evaluator.add_plain_inplace(o1_3_0Cip[j], half_3Plain);
			evaluator.multiply_inplace(layer_1_3_deltaCip[0][j], o1_3_0Cip[j]);
			evaluator.relinearize_inplace(layer_1_3_deltaCip[0][j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_3_deltaCip[0][j]);
		}
		//研究者：将8个密文的后半部分广播
		//数据提供者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < 8 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[0][j - 1][k][m] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
				}
				ckks_encoder.encode(backnoises[0][j - 1][k], levelFour, levelFourScale, tempPlain);
				decryptor[j]->partial_decrypt(layer_1_3_deltaCip[0][k], tempPlain, pDecPlain[k][j - 1]);
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[0][j - 1][k][m] *= -27.902113901977722956 / pow(2, 15.5);
				}
				ckks_encoder.encode(backnoises[0][j - 1][k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesOriCip[j - 1][k]);
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[0][j - 1][k][m] /= 5.282245157314995551;
				}
			}
			dataCombineEncode(temps, backnoises[0][j - 1], backnoises[1][j - 1], batchsize * time / 8, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), empty, 16 * (int)1.000);
			for (int k = 0; k < 16 * (int)1.000; k++)
			{
				ckks_encoder.encode(temps[k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesCip[j - 1][k]);
			}
		}
		//数据提供者：将(16 + 8)×(partyNum - 1)个噪声密文发送，8×(partyNum - 1)个部分解密发送
		//研究者：完全解密并重新编码(将8个layer_1_3_deltaOriCip[0]发送给数据提供者)
		for (int j = 0; j < 8 * (int)1.000; j++)
		{
			decryptor[0]->final_decrypt(layer_1_3_deltaCip[0][j], pDecPlain[j], partyNum - 1, tempPlain);
			ckks_encoder.decode(tempPlain, noises[0][j]);
		}
		for (int j = 0; j < 8 * (int)1.000; j++)
		{
			for (int k = 0; k < slotCount; k++)
			{
				noises[0][j][k] *= 27.902113901977722956 / pow(2, 15.5);
			}
			ckks_encoder.encode(noises[0][j], levelTwo, levelTwoScale, tempPlain);
			layer_1_3_deltaOriCip[0][j][0] = noisesOriCip[0][j];
			evaluator.add_plain_inplace(layer_1_3_deltaOriCip[0][j][0], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_3_deltaOriCip[0][j][0], noisesOriCip[k - 1][j]);
				k++;
			}
			for (int k = 0; k < slotCount; k++)
			{
				noises[0][j][k] /= 5.282245157314995551;
			}
		}
		dataCombineEncode(temps, noises[0], noises[1], batchsize * time / 8, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), empty, 16 * (int)1.000);
		for (int j = 0; j < 16 * (int)1.000; j++)
		{
			ckks_encoder.encode(temps[j], levelTwo, levelTwoScale, tempPlain);
			evaluator.add_plain_inplace(noisesCip[0][j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(noisesCip[0][j], noisesCip[k - 1][j]);
				k++;
			}
		}
		//研究者：将16个密文发送
		//数据提供者：反向传播8倍数_2
		for (int j = 0; j < 4; j++)
		{
			evaluator.multiply_plain(wb_3TranCip, whb_filterPlain[j], wb_3EncodeCip[j]);
			for (int k = pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(wb_3EncodeCip[j], -k, galois_keys, tempCip);
				evaluator.add_inplace(wb_3EncodeCip[j], tempCip);
			}
			if (j != 0)
			{
				evaluator.rotate_vector_inplace(wb_3EncodeCip[j], j * pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4, galois_keys);
			}
			evaluator.rescale_to_next_inplace(wb_3EncodeCip[j]);
		}
		parameterRotate(wb_3EncodeCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), paramCip[0], evaluator, galois_keys, one_1Plain);
		matrixMultiple_CipCip(paramCip[0], noisesCip[0], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 16 * (int)1.000, layer_1_2_deltaCip[1], evaluator, galois_keys);
		for (int j = 0; j < 16 * (int)1.000; j++)
		{
			evaluator.relinearize_inplace(layer_1_2_deltaCip[1][j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_2_deltaCip[1][j]);
			evaluator.add_plain_inplace(o1_2_1Cip[j], half_3Plain);
			evaluator.multiply_inplace(layer_1_2_deltaCip[1][j], o1_2_1Cip[j]);
			evaluator.relinearize_inplace(layer_1_2_deltaCip[1][j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_2_deltaCip[1][j]);
		}
		//数据提供者：将16个密文的后半部分广播
		//研究者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < 16 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[1][j - 1][k][m] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
				}
				ckks_encoder.encode(backnoises[1][j - 1][k], levelFour, levelFourScale, tempPlain);
				decryptor[j]->partial_decrypt(layer_1_2_deltaCip[1][k], tempPlain, pDecPlain[k][j - 1]);
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[1][j - 1][k][m] *= -27.902113901977722956 / pow(2, 15.5);
				}
				ckks_encoder.encode(backnoises[1][j - 1][k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesOriCip[j - 1][k]);
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[1][j - 1][k][m] /= 5.282245157314995551;
				}
			}
			matrixRefresh(backnoises[1][j - 1], temps, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 16 * (int)1.000, batchsize * time / 8);
			for (int k = 0; k < 16 * (int)1.000; k++)
			{
				ckks_encoder.encode(temps[k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesCip[j - 1][k]);
			}
		}
		//研究者：将16×2×(partyNum - 1)个噪声密文发送，16×(partyNum - 1)个部分解密发送
		//数据提供者：完全解密并重新编码
		for (int j = 0; j < 16 * (int)1.000; j++)
		{
			decryptor[0]->final_decrypt(layer_1_2_deltaCip[1][j], pDecPlain[j], partyNum - 1, tempPlain);
			ckks_encoder.decode(tempPlain, noises[1][j]);
		}
		for (int j = 0; j < 16 * (int)1.000; j++)
		{
			for (int k = 0; k < slotCount; k++)
			{
				noises[1][j][k] *= 27.902113901977722956 / pow(2, 15.5);
			}
			ckks_encoder.encode(noises[1][j], levelTwo, levelTwoScale, tempPlain);
			layer_1_2_deltaOriCip[1][j][0] = noisesOriCip[0][j];
			evaluator.add_plain_inplace(layer_1_2_deltaOriCip[1][j][0], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_2_deltaOriCip[1][j][0], noisesOriCip[k - 1][j]);
				k++;
			}
			for (int k = 0; k < slotCount; k++)
			{
				noises[1][j][k] /= 5.282245157314995551;
			}
		}
		matrixRefresh(noises[1], temps, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 16 * (int)1.000, batchsize * time / 8);
		for (int j = 0; j < 16 * (int)1.000; j++)
		{
			ckks_encoder.encode(temps[j], levelTwo, levelTwoScale, tempPlain);
			layer_1_2_deltaCip[1][j] = noisesCip[0][j];
			evaluator.add_plain_inplace(layer_1_2_deltaCip[1][j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_2_deltaCip[1][j], noisesCip[k - 1][j]);
				k++;
			}
		}
		//数据提供者：反向传播8余4_2
		for (int j = 0; j < 4; j++)
		{
			evaluator.multiply_plain(wh_2TranCip, whb_filterPlain[j], wh_2EncodeCip[j]);
			for (int k = pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(wh_2EncodeCip[j], -k, galois_keys, tempCip);
				evaluator.add_inplace(wh_2EncodeCip[j], tempCip);
			}
			if (j != 0)
			{
				evaluator.rotate_vector_inplace(wh_2EncodeCip[j], j * pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4, galois_keys);
			}
			evaluator.rescale_to_next_inplace(wh_2EncodeCip[j]);
		}
		parameterRotate(wh_2EncodeCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), paramCip[1], evaluator, galois_keys, one_1Plain);
		matrixMultiple_CipCip(paramCip[1], layer_1_2_deltaCip[1], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 16 * (int)1.000, layer_1_2_deltaCip[0], evaluator, galois_keys);
		for (int j = 0; j < 16 * (int)1.000; j++)
		{
			evaluator.relinearize_inplace(layer_1_2_deltaCip[0][j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_2_deltaCip[0][j]);
			evaluator.add_plain_inplace(o1_2_0Cip[j], half_3Plain);
			evaluator.multiply_inplace(layer_1_2_deltaCip[0][j], o1_2_0Cip[j]);
			evaluator.relinearize_inplace(layer_1_2_deltaCip[0][j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_2_deltaCip[0][j]);
		}
		//数据提供者：将16个密文的后半部分广播
		//研究者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < 16 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[0][j - 1][k][m] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
				}
				ckks_encoder.encode(backnoises[0][j - 1][k], levelFour, levelFourScale, tempPlain);
				decryptor[j]->partial_decrypt(layer_1_2_deltaCip[0][k], tempPlain, pDecPlain[k][j - 1]);
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[0][j - 1][k][m] *= -27.902113901977722956 / pow(2, 15.5);
				}
				ckks_encoder.encode(backnoises[0][j - 1][k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesOriCip[j - 1][k]);
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[0][j - 1][k][m] /= 5.282245157314995551;
				}
			}
			dataCombineEncode(temps, backnoises[0][j - 1], backnoises[1][j - 1], batchsize * time / 4, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), empty, 28 * (int)1.000);
			for (int k = 0; k < 28 * (int)1.000; k++)
			{
				ckks_encoder.encode(temps[k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesCip[j - 1][k]);
			}
		}
		//研究者：将(28 + 16)×(partyNum - 1)个噪声密文发送，16×(partyNum - 1)个部分解密发送
		//数据提供者：完全解密并重新编码
		for (int j = 0; j < 16 * (int)1.000; j++)
		{
			decryptor[0]->final_decrypt(layer_1_2_deltaCip[0][j], pDecPlain[j], partyNum - 1, tempPlain);
			ckks_encoder.decode(tempPlain, noises[0][j]);
		}
		for (int j = 0; j < 16 * (int)1.000; j++)
		{
			for (int k = 0; k < slotCount; k++)
			{
				noises[0][j][k] *= 27.902113901977722956 / pow(2, 15.5);
			}
			ckks_encoder.encode(noises[0][j], levelTwo, levelTwoScale, tempPlain);
			layer_1_2_deltaOriCip[0][j][0] = noisesOriCip[0][j];
			evaluator.add_plain_inplace(layer_1_2_deltaOriCip[0][j][0], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_2_deltaOriCip[0][j][0], noisesOriCip[k - 1][j]);
				k++;
			}
			for (int k = 0; k < slotCount; k++)
			{
				noises[0][j][k] /= 5.282245157314995551;
			}
		}
		dataCombineEncode(temps, noises[0], noises[1], batchsize * time / 4, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), empty, 28 * (int)1.000);
		for (int j = 0; j < 28 * (int)1.000; j++)
		{
			ckks_encoder.encode(temps[j], levelTwo, levelTwoScale, tempPlain);
			evaluator.add_plain_inplace(noisesCip[0][j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(noisesCip[0][j], noisesCip[k - 1][j]);
				k++;
			}
		}
		//数据提供者：将28个密文发送
		//研究者：反向传播4倍数_1
		for (int j = 0; j < 4; j++)
		{
			evaluator.multiply_plain(wb_2TranCip, whb_filterPlain[j], wb_2EncodeCip[j]);
			for (int k = pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(wb_2EncodeCip[j], -k, galois_keys, tempCip);
				evaluator.add_inplace(wb_2EncodeCip[j], tempCip);
			}
			if (j != 0)
			{
				evaluator.rotate_vector_inplace(wb_2EncodeCip[j], j * pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4, galois_keys);
			}
			evaluator.rescale_to_next_inplace(wb_2EncodeCip[j]);
		}
		parameterRotate(wb_2EncodeCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), paramCip[0], evaluator, galois_keys, one_1Plain);
		matrixMultiple_CipCip(paramCip[0], noisesCip[0], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 28 * (int)1.000, layer_1_1_deltaCip[1], evaluator, galois_keys);
		for (int j = 0; j < 28 * (int)1.000; j++)
		{
			evaluator.relinearize_inplace(layer_1_1_deltaCip[1][j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_1_deltaCip[1][j]);
			evaluator.add_plain_inplace(o1_1_1Cip[j], half_3Plain);
			evaluator.multiply_inplace(layer_1_1_deltaCip[1][j], o1_1_1Cip[j]);
			evaluator.relinearize_inplace(layer_1_1_deltaCip[1][j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_1_deltaCip[1][j]);
		}
		//研究者：将28个密文的后半部分广播
		//数据提供者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < 28 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[1][j - 1][k][m] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
				}
				ckks_encoder.encode(backnoises[1][j - 1][k], levelFour, levelFourScale, tempPlain);
				decryptor[j]->partial_decrypt(layer_1_1_deltaCip[1][k], tempPlain, pDecPlain[k][j - 1]);
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[1][j - 1][k][m] *= -27.902113901977722956 / pow(2, 15.5);
				}
				ckks_encoder.encode(backnoises[1][j - 1][k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesOriCip[j - 1][k]);
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[1][j - 1][k][m] /= 5.282245157314995551;
				}
			}
			matrixRefresh(backnoises[1][j - 1], temps, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 28 * (int)1.000, batchsize * time / 4);
			for (int k = 0; k < 28 * (int)1.000; k++)
			{
				ckks_encoder.encode(temps[k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesCip[j - 1][k]);
			}
		}
		//数据提供者：将28×2×(partyNum - 1)个噪声密文发送，28×(partyNum - 1)个部分解密发送
		//研究者：完全解密并重新编码(将28个layer_1_1_deltaOriCip[1]发送给数据提供者)
		for (int j = 0; j < 28 * (int)1.000; j++)
		{
			decryptor[0]->final_decrypt(layer_1_1_deltaCip[1][j], pDecPlain[j], partyNum - 1, tempPlain);
			ckks_encoder.decode(tempPlain, noises[1][j]);
		}
		for (int j = 0; j < 28 * (int)1.000; j++)
		{
			for (int k = 0; k < slotCount; k++)
			{
				noises[1][j][k] *= 27.902113901977722956 / pow(2, 15.5);
			}
			ckks_encoder.encode(noises[1][j], levelTwo, levelTwoScale, tempPlain);
			layer_1_1_deltaOriCip[1][j][0] = noisesOriCip[0][j];
			evaluator.add_plain_inplace(layer_1_1_deltaOriCip[1][j][0], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_1_deltaOriCip[1][j][0], noisesOriCip[k - 1][j]);
				k++;
			}
			for (int k = 0; k < slotCount; k++)
			{
				noises[1][j][k] /= 5.282245157314995551;
			}
		}
		matrixRefresh(noises[1], temps, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 28 * (int)1.000, batchsize * time / 4);
		for (int j = 0; j < 28 * (int)1.000; j++)
		{
			ckks_encoder.encode(temps[j], levelTwo, levelTwoScale, tempPlain);
			layer_1_1_deltaCip[1][j] = noisesCip[0][j];
			evaluator.add_plain_inplace(layer_1_1_deltaCip[1][j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_1_deltaCip[1][j], noisesCip[k - 1][j]);
				k++;
			}
		}
		//研究者：反向传播4余2_1
		for (int j = 0; j < 4; j++)
		{
			evaluator.multiply_plain(wh_1TranCip, whb_filterPlain[j], wh_1EncodeCip[j]);
			for (int k = pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(wh_1EncodeCip[j], -k, galois_keys, tempCip);
				evaluator.add_inplace(wh_1EncodeCip[j], tempCip);
			}
			if (j != 0)
			{
				evaluator.rotate_vector_inplace(wh_1EncodeCip[j], j * pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4, galois_keys);
			}
			evaluator.rescale_to_next_inplace(wh_1EncodeCip[j]);
		}
		parameterRotate(wh_1EncodeCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), paramCip[1], evaluator, galois_keys, one_1Plain);
		matrixMultiple_CipCip(paramCip[1], layer_1_1_deltaCip[1], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 28 * (int)1.000, layer_1_1_deltaCip[0], evaluator, galois_keys);
		for (int j = 0; j < 28 * (int)1.000; j++)
		{
			evaluator.relinearize_inplace(layer_1_1_deltaCip[0][j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_1_deltaCip[0][j]);
			evaluator.add_plain_inplace(o1_1_0Cip[j], half_3Plain);
			evaluator.multiply_inplace(layer_1_1_deltaCip[0][j], o1_1_0Cip[j]);
			evaluator.relinearize_inplace(layer_1_1_deltaCip[0][j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_1_deltaCip[0][j]);
		}
		//研究者：将28个密文的后半部分广播
		//数据提供者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < 28 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[0][j - 1][k][m] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
				}
				ckks_encoder.encode(backnoises[0][j - 1][k], levelFour, levelFourScale, tempPlain);
				decryptor[j]->partial_decrypt(layer_1_1_deltaCip[0][k], tempPlain, pDecPlain[k][j - 1]);
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[0][j - 1][k][m] *= -27.902113901977722956 / pow(2, 15.5);
				}
				ckks_encoder.encode(backnoises[0][j - 1][k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesOriCip[j - 1][k]);
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[0][j - 1][k][m] /= 5.282245157314995551;
				}
			}
			dataCombineEncode(temps, backnoises[0][j - 1], backnoises[1][j - 1], batchsize * time / 2, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), empty, 52 * (int)1.000);
			for (int k = 0; k < 52 * (int)1.000; k++)
			{
				ckks_encoder.encode(temps[k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesCip[j - 1][k]);
			}
		}
		//数据提供者：将(52 + 28)×(partyNum - 1)个噪声密文发送，28×(partyNum - 1)个部分解密发送
		//研究者：完全解密并重新编码(将28个layer_1_1_deltaOriCip[1]发送给数据提供者)
		for (int j = 0; j < 28 * (int)1.000; j++)
		{
			decryptor[0]->final_decrypt(layer_1_1_deltaCip[0][j], pDecPlain[j], partyNum - 1, tempPlain);
			ckks_encoder.decode(tempPlain, noises[0][j]);
		}
		for (int j = 0; j < 28 * (int)1.000; j++)
		{
			for (int k = 0; k < slotCount; k++)
			{
				noises[0][j][k] *= 27.902113901977722956 / pow(2, 15.5);
			}
			ckks_encoder.encode(noises[0][j], levelTwo, levelTwoScale, tempPlain);
			layer_1_1_deltaOriCip[0][j][0] = noisesOriCip[0][j];
			evaluator.add_plain_inplace(layer_1_1_deltaOriCip[0][j][0], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_1_deltaOriCip[0][j][0], noisesOriCip[k - 1][j]);
				k++;
			}
			for (int k = 0; k < slotCount; k++)
			{
				noises[0][j][k] /= 5.282245157314995551;
			}
		}
		dataCombineEncode(temps, noises[0], noises[1], batchsize * time / 2, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), empty, 52 * (int)1.000);
		for (int j = 0; j < 52 * (int)1.000; j++)
		{
			ckks_encoder.encode(temps[j], levelTwo, levelTwoScale, tempPlain);
			evaluator.add_plain_inplace(noisesCip[0][j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(noisesCip[0][j], noisesCip[k - 1][j]);
				k++;
			}
		}
		//研究者：将52个密文发送
		//数据提供者：反向传播2倍数_0
		for (int j = 0; j < 4; j++)
		{
			evaluator.multiply_plain(wb_1TranCip, whb_filterPlain[j], wb_1EncodeCip[j]);
			for (int k = pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(wb_1EncodeCip[j], -k, galois_keys, tempCip);
				evaluator.add_inplace(wb_1EncodeCip[j], tempCip);
			}
			if (j != 0)
			{
				evaluator.rotate_vector_inplace(wb_1EncodeCip[j], j * pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4, galois_keys);
			}
			evaluator.rescale_to_next_inplace(wb_1EncodeCip[j]);
		}
		parameterRotate(wb_1EncodeCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), paramCip[0], evaluator, galois_keys, one_1Plain);
		matrixMultiple_CipCip(paramCip[0], noisesCip[0], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 52 * (int)1.000, layer_1_0_deltaCip[1], evaluator, galois_keys);
		for (int j = 0; j < 52 * (int)1.000; j++)
		{
			evaluator.relinearize_inplace(layer_1_0_deltaCip[1][j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_0_deltaCip[1][j]);
			evaluator.add_plain_inplace(o1_0_1Cip[j], half_3Plain);
			evaluator.multiply_inplace(layer_1_0_deltaCip[1][j], o1_0_1Cip[j]);
			evaluator.relinearize_inplace(layer_1_0_deltaCip[1][j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_0_deltaCip[1][j]);
		}
		//数据提供者：将52个密文的后半部分广播
		//研究者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < 52 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[1][j - 1][k][m] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
				}
				ckks_encoder.encode(backnoises[1][j - 1][k], levelFour, levelFourScale, tempPlain);
				decryptor[j]->partial_decrypt(layer_1_0_deltaCip[1][k], tempPlain, pDecPlain[k][j - 1]);
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[1][j - 1][k][m] *= -27.902113901977722956 / pow(2, 15.5);
				}
				ckks_encoder.encode(backnoises[1][j - 1][k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesOriCip[j - 1][k]);
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[1][j - 1][k][m] /= 5.282245157314995551;
				}
			}
			matrixRefresh(backnoises[1][j - 1], temps, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 52 * (int)1.000, batchsize * time / 2);
			for (int k = 0; k < 52 * (int)1.000; k++)
			{
				ckks_encoder.encode(temps[k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesCip[j - 1][k]);
			}
		}
		//研究者：将52×2×(partyNum - 1)个噪声密文发送，52×(partyNum - 1)个部分解密发送
		//数据提供者：完全解密并重新编码
		for (int j = 0; j < 52 * (int)1.000; j++)
		{
			decryptor[0]->final_decrypt(layer_1_0_deltaCip[1][j], pDecPlain[j], partyNum - 1, tempPlain);
			ckks_encoder.decode(tempPlain, noises[1][j]);
		}
		for (int j = 0; j < 52 * (int)1.000; j++)
		{
			for (int k = 0; k < slotCount; k++)
			{
				noises[1][j][k] *= 27.902113901977722956 / pow(2, 15.5);
			}
			ckks_encoder.encode(noises[1][j], levelTwo, levelTwoScale, tempPlain);
			layer_1_0_deltaOriCip[1][j][0] = noisesOriCip[0][j];
			evaluator.add_plain_inplace(layer_1_0_deltaOriCip[1][j][0], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_0_deltaOriCip[1][j][0], noisesOriCip[k - 1][j]);
				k++;
			}
			for (int k = 0; k < slotCount; k++)
			{
				noises[1][j][k] /= 5.282245157314995551;
			}
		}
		matrixRefresh(noises[1], temps, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 52 * (int)1.000, batchsize * time / 2);
		for (int j = 0; j < 52 * (int)1.000; j++)
		{
			ckks_encoder.encode(temps[j], levelTwo, levelTwoScale, tempPlain);
			layer_1_0_deltaCip[1][j] = noisesCip[0][j];
			evaluator.add_plain_inplace(layer_1_0_deltaCip[1][j], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_0_deltaCip[1][j], noisesCip[k - 1][j]);
				k++;
			}
		}
		//数据提供者：反向传播2余1_0
		for (int j = 0; j < 4; j++)
		{
			evaluator.multiply_plain(wh_0TranCip, whb_filterPlain[j], wh_0EncodeCip[j]);
			for (int k = pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(wh_0EncodeCip[j], -k, galois_keys, tempCip);
				evaluator.add_inplace(wh_0EncodeCip[j], tempCip);
			}
			if (j != 0)
			{
				evaluator.rotate_vector_inplace(wh_0EncodeCip[j], j * pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4, galois_keys);
			}
			evaluator.rescale_to_next_inplace(wh_0EncodeCip[j]);
		}
		parameterRotate(wh_0EncodeCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), paramCip[1], evaluator, galois_keys, one_1Plain);
		matrixMultiple_CipCip(paramCip[1], layer_1_0_deltaCip[1], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), 52 * (int)1.000, layer_1_0_deltaCip[0], evaluator, galois_keys);
		for (int j = 0; j < 52 * (int)1.000; j++)
		{
			evaluator.relinearize_inplace(layer_1_0_deltaCip[0][j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_0_deltaCip[0][j]);
			evaluator.add_plain_inplace(o1_0_0Cip[j], half_3Plain);
			evaluator.multiply_inplace(layer_1_0_deltaCip[0][j], o1_0_0Cip[j]);
			evaluator.relinearize_inplace(layer_1_0_deltaCip[0][j], relin_keys);
			evaluator.rescale_to_next_inplace(layer_1_0_deltaCip[0][j]);
		}
		//数据提供者：将52个密文的后半部分广播
		//研究者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < 52 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[0][j - 1][k][m] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
				}
				ckks_encoder.encode(backnoises[0][j - 1][k], levelFour, levelFourScale, tempPlain);
				decryptor[j]->partial_decrypt(layer_1_0_deltaCip[0][k], tempPlain, pDecPlain[k][j - 1]);
				for (int m = 0; m < slotCount; m++)
				{
					backnoises[0][j - 1][k][m] *= -27.902113901977722956 / pow(2, 15.5);
				}
				ckks_encoder.encode(backnoises[0][j - 1][k], levelTwo, levelTwoScale, tempPlain);
				encryptor.encrypt(tempPlain, noisesOriCip[j - 1][k]);
			}
		}
		//研究者：将52×(partyNum - 1)个噪声密文发送，52×(partyNum - 1)个部分解密发送
		//数据提供者：完全解密并重新编码
		for (int j = 0; j < 52 * (int)1.000; j++)
		{
			decryptor[0]->final_decrypt(layer_1_0_deltaCip[0][j], pDecPlain[j], partyNum - 1, tempPlain);
			ckks_encoder.decode(tempPlain, noises[0][j]);
		}
		for (int j = 0; j < 52 * (int)1.000; j++)
		{
			for (int k = 0; k < slotCount; k++)
			{
				noises[0][j][k] *= 27.902113901977722956 / pow(2, 15.5);
			}
			ckks_encoder.encode(noises[0][j], levelTwo, levelTwoScale, tempPlain);
			layer_1_0_deltaOriCip[0][j][0] = noisesOriCip[0][j];
			evaluator.add_plain_inplace(layer_1_0_deltaOriCip[0][j][0], tempPlain);
			int k = 2;
			while (k < partyNum)
			{
				evaluator.add_inplace(layer_1_0_deltaOriCip[0][j][0], noisesOriCip[k - 1][j]);
				k++;
			}
		}
		//数据提供者：保存隐含层和输出层之间的梯度
		for (int j = 0; j < 26 * (int)1.000; j++)
		{
			coeff[j] = empty;
		}
		Pos = 0;
		for (int j = 0; j < batchsize; j++)
		{
			for (int k = 0; k < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k++)
			{
				coeff[Pos + j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) % 2][j % (int)((slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) - pow(2, ceil(log2(hidenode * hidenode / 4) / 2))) * 2) / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2 + 1)) * pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j % (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + (long long)k * slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2))] = (outcome[(i * batchsize + j) % row] * 0.84 + 0.16) * alpha;
			}
			if (j % (int)((slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) - pow(2, ceil(log2(hidenode * hidenode / 4) / 2))) * 2) == (slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) - pow(2, ceil(log2(hidenode * hidenode / 4) / 2))) * 2 - 1)
			{
				Pos += 2;
			}
		}
		for (int j = 0; j < 2 * (int)1.000; j++)
		{
			ckks_encoder.encode(coeff[j], levelTwo, levelTwoScale, tempsPlain[j]);
		}
		for (int j = 0; j < 2 * (int)1.000; j++)
		{
			evaluator.multiply_plain(layer_2_deltaCip[j], tempsPlain[j], tempCip);
			evaluator.rescale_to_next_inplace(tempCip);
			for (int k = 0; k < 2; k++)
			{
				evaluator.multiply(layer_1_4_2Cip[k * 2 + j + j / 2 * 2], tempCip, tempsCip[k * 2 + j + j / 2 * 2]);
				evaluator.relinearize_inplace(tempsCip[k * 2 + j + j / 2 * 2], relin_keys);
			}
		}
		for (int j = 0; j < 2 * (int)1.000; j++)
		{
			evaluator.add_inplace(tempsCip[j * 2], tempsCip[j * 2 + 1]);
			for (int k = 1; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(tempsCip[j * 2], k, galois_keys, tempCip);
				evaluator.add_inplace(tempsCip[j * 2], tempCip);
			}
			evaluator.rescale_to_next_inplace(tempsCip[j * 2]);
		}
		//数据提供者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < hidenode; k++)
			{
				w1_delta[k] = 0;
			}
			b1_delta = 0;
			for (int k = 0; k < 2 * (int)1.000; k++)
			{
				for (int m = 0; m < slotCount; m++)
				{
					noises[j - 1][k][m] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
				}
				ckks_encoder.encode(noises[j - 1][k], levelFour, levelFourScale, tempPlain);
				decryptor[j]->partial_decrypt(tempsCip[k * 2], tempPlain, pDecPlain[k][j - 1]);
			}
			for (int k = 0; k < 2 * (int)1.000; k++)
			{
				for (int m = 0; m < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); m++)
				{
					if (k * pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m < hidenode)
					{
						w1_delta[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m] -= noises[j - 1][k][(long long)m * slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2))] / pow(2, 23.25 - 15.5 / 6.0);
					}
					else
					{
						b1_delta -= noises[j - 1][k][(long long)m * slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2))] / pow(2, 23.25 - 15.5 / 6.0);
					}
				}
			}
			noisesCip[j - 1][0] = zeroCip;
			outParameterEncode(noisesCip[j - 1][0], hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w1_delta, b1_delta);
			tempCip = noisesCip[j - 1][0];
			evaluator.add_inplace(noisesCip[j - 1][0], wb1_NoiseCip[j - 1]);
			evaluator.add_inplace(wb1_NoiseCip[j - 1], tempCip);
		}
		//数据提供者：将partyNum - 1个噪声密文发送，2×(partyNum - 1)个部分解密发送
		//研究者：完全解密并重新编码
		for (int j = 0; j < 2 * (int)1.000; j++)
		{
			decryptor[0]->final_decrypt(tempsCip[j * 2], pDecPlain[j], partyNum - 1, tempPlain);
			ckks_encoder.decode(tempPlain, noises[1][j]);
		}
		for (int j = 0; j < 2 * (int)1.000; j++)
		{
			for (int k = 0; k < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k++)
			{
				if (j * pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k < hidenode)
				{
					w1[j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k] += noises[1][j][(long long)k * slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2))] / pow(2, 23.25 - 15.5 / 6.0);
				}
				else
				{
					b1 += noises[1][j][(long long)k * slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2))] / pow(2, 23.25 - 15.5 / 6.0);
				}
			}
		}
		wb1_Cip = zeroCip;
		outParameterEncode(wb1_Cip, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w1, b1);
		for (int j = 1; j < partyNum; j++)
		{
			evaluator.add_inplace(wb1_Cip, noisesCip[j - 1][0]);
		}
		//数据提供者：保存隐含层之间的梯度(48_4)
		for (int j = 0; j < 2 * (int)1.000; j++)
		{
			ckks_encoder.encode(coeff[j], levelOne, levelOneScale, tempsPlain[j]);
		}
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 4 * (int)1.000; k++)
			{
				for (int m = 1; m < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); m++)
				{
					evaluator.rotate_vector(layer_1_4_deltaOriCip[j][k][m - 1], -slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), galois_keys, layer_1_4_deltaOriCip[j][k][m]);
				}
			}
		}
		for (int j = 0; j < 1 * (int)1.000; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				evaluator.multiply_plain(layer_1_4_1OriCip[k + j * 4], tempsPlain[k % 2 + j * 2], tempCip);
				evaluator.rescale_to_next_inplace(tempCip);
				for (int m = 0; m < 2; m++)
				{
					for (int n = 0; n < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); n++)
					{
						evaluator.multiply(layer_1_4_deltaOriCip[2][m * 2 + k % 2 + j * 4][n], tempCip, tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
						evaluator.relinearize_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n], relin_keys);
						evaluator.rescale_to_next_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
					}
				}
			}
		}
		for (int j = 0; j < 1 * (int)1.000; j++)
		{
			for (int k = 0; k < 64; k++)
			{
				evaluator.add_inplace(tempsCip[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8], tempsCip[(k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + 1) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8]);
			}
		}
		for (int j = 0; j < (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 * (int)1.000; j++)
		{
			tempsCip[j + (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 * (int)1.000] = tempsCip[j];
		}
		//数据提供者：保存隐含层之间的梯度(32_4)
		for (int j = 0; j < 1 * (int)1.000; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				evaluator.multiply_plain(layer_1_4_0OriCip[k + j * 4], tempsPlain[k % 2 + j * 2], tempCip);
				evaluator.rescale_to_next_inplace(tempCip);
				for (int m = 0; m < 2; m++)
				{
					for (int n = 0; n < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); n++)
					{
						evaluator.multiply(layer_1_4_deltaOriCip[1][m * 2 + k % 2 + j * 4][n], tempCip, tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
						evaluator.relinearize_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n], relin_keys);
						evaluator.rescale_to_next_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
					}
				}
			}
		}
		for (int j = 0; j < 1 * (int)1.000; j++)
		{
			for (int k = 0; k < 64; k++)
			{
				evaluator.add_inplace(tempsCip[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8], tempsCip[(k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + 1) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8]);
			}
		}
		for (int j = 0; j < (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 * (int)1.000; j++)
		{
			evaluator.add_inplace(tempsCip[j], tempsCip[j + (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 * (int)1.000]);
		}
		for (int j = 0; j < 64; j++)
		{
			for (int k = 1; k < 1 * (int)1.000; k++)
			{
				evaluator.add_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j + (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 * k]);
			}
			for (int k = 1; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], k, galois_keys, tempCip);
				evaluator.add_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], tempCip);
			}
			evaluator.multiply_plain_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], one_filterPlain);
		}
		for (int j = 1; j < 64; j++)
		{
			evaluator.rotate_vector_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], -updatestep[j], galois_keys);
			evaluator.add_inplace(tempsCip[0], tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j]);
		}
		wh_4_deltaCip = tempsCip[0];
		evaluator.rescale_to_next_inplace(wh_4_deltaCip);
		//数据提供者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < hidenode; k++)
			{
				for (int m = 0; m < hidenode; m++)
				{
					wh_4_delta[k][m] = 0;
				}
			}
			for (int k = 0; k < slotCount; k++)
			{
				noises[j - 1][0][k] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
			}
			ckks_encoder.encode(noises[j - 1][0], levelFour, levelFourScale, tempPlain);
			decryptor[j]->partial_decrypt(wh_4_deltaCip, tempPlain, pDecPlain[0][j - 1]);
			for (int k = 0; k < hidenode; k++)
			{
				for (int m = 0; m < hidenode; m++)
				{
					wh_4_delta[m][k] -= noises[j - 1][0][(long long)m * 256 % 4096 + updatestep[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 2 + m / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(m - k, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
				}
			}
			noisesCip[j - 1][0] = zeroCip;
			parameterEncode(noisesCip[j - 1][0], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_4_delta);
			tempCip = noisesCip[j - 1][0];
			evaluator.add_inplace(noisesCip[j - 1][0], wh_4_NoiseCip[j - 1]);
			evaluator.add_inplace(wh_4_NoiseCip[j - 1], tempCip);
			noisesCip[j - 1][1] = zeroCip;
			parameterEncodeTranspose(noisesCip[j - 1][1], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_4_delta);
			tempCip = noisesCip[j - 1][1];
			evaluator.add_inplace(noisesCip[j - 1][1], wh_4_TranNoiseCip[j - 1]);
			evaluator.add_inplace(wh_4_TranNoiseCip[j - 1], tempCip);
		}
		//数据提供者：将2×(partyNum - 1)个噪声密文发送，partyNum - 1个部分解密发送
		//研究者：完全解密并重新编码
		decryptor[0]->final_decrypt(wh_4_deltaCip, pDecPlain[0], partyNum - 1, tempPlain);
		ckks_encoder.decode(tempPlain, noises[1][0]);
		for (int j = 0; j < hidenode; j++)
		{
			for (int k = 0; k < hidenode; k++)
			{
				wh_4[k][j] += noises[1][0][(long long)k * 256 % 4096 + updatestep[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 2 + k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(k - j, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
			}
		}
		wh_4Cip = zeroCip;
		parameterEncode(wh_4Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_4);
		for (int j = 1; j < partyNum; j++)
		{
			evaluator.add_inplace(wh_4Cip, noisesCip[j - 1][0]);
		}
		wh_4TranCip = zeroCip;
		parameterEncodeTranspose(wh_4TranCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_4);
		for (int j = 1; j < partyNum; j++)
		{
			evaluator.add_inplace(wh_4TranCip, noisesCip[j - 1][1]);
		}
		//数据提供者：保存输入层和隐含层之间的梯度(48_4)
		for (int j = 0; j < 1 * (int)1.000; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				evaluator.multiply_plain(layer_1_3_1OriCip[k + j * 4 + 8 * (int)1.000], tempsPlain[k % 2 + j * 2], tempCip);
				evaluator.rescale_to_next_inplace(tempCip);
				for (int m = 0; m < 2; m++)
				{
					for (int n = 0; n < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); n++)
					{
						evaluator.multiply(layer_1_4_deltaOriCip[2][m * 2 + k % 2 + j * 4][n], tempCip, tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
						evaluator.relinearize_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n], relin_keys);
						evaluator.rescale_to_next_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
					}
				}
			}
		}
		for (int j = 0; j < 1 * (int)1.000; j++)
		{
			for (int k = 0; k < 64; k++)
			{
				evaluator.add_inplace(tempsCip[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8], tempsCip[(k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + 1) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8]);
			}
		}
		for (int j = 0; j < (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 * (int)1.000; j++)
		{
			tempsCip[j + (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 * (int)1.000 * 2] = tempsCip[j];
		}
		//数据提供者：保存输入层和隐含层之间的梯度(32_4)
		for (int j = 0; j < 1 * (int)1.000; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				evaluator.multiply_plain(layer_1_3_1OriCip[k + j * 4 + 4 * (int)1.000], tempsPlain[k % 2 + j * 2], tempCip);
				evaluator.rescale_to_next_inplace(tempCip);
				for (int m = 0; m < 2; m++)
				{
					for (int n = 0; n < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); n++)
					{
						evaluator.multiply(layer_1_4_deltaOriCip[1][m * 2 + k % 2 + j * 4][n], tempCip, tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
						evaluator.relinearize_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n], relin_keys);
						evaluator.rescale_to_next_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
					}
				}
			}
		}
		for (int j = 0; j < 1 * (int)1.000; j++)
		{
			for (int k = 0; k < 64; k++)
			{
				evaluator.add_inplace(tempsCip[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8], tempsCip[(k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + 1) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8]);
			}
		}
		for (int j = 0; j < (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 * (int)1.000; j++)
		{
			tempsCip[j + (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 * (int)1.000] = tempsCip[j];
		}
		//数据提供者：保存输入层和隐含层之间的梯度(16_4)
		for (int j = 0; j < 1 * (int)1.000; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				evaluator.multiply_plain(layer_1_3_1OriCip[k + j * 4], tempsPlain[k % 2 + j * 2], tempCip);
				evaluator.rescale_to_next_inplace(tempCip);
				for (int m = 0; m < 2; m++)
				{
					for (int n = 0; n < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); n++)
					{
						evaluator.multiply(layer_1_4_deltaOriCip[0][m * 2 + k % 2 + j * 4][n], tempCip, tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
						evaluator.relinearize_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n], relin_keys);
						evaluator.rescale_to_next_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
					}
				}
			}
		}
		for (int j = 0; j < 1 * (int)1.000; j++)
		{
			for (int k = 0; k < 64; k++)
			{
				evaluator.add_inplace(tempsCip[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8], tempsCip[(k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + 1) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8]);
			}
		}
		for (int j = 0; j < (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 * (int)1.000; j++)
		{
			evaluator.add_inplace(tempsCip[j], tempsCip[j + (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 * (int)1.000 * 2]);
			evaluator.add_inplace(tempsCip[j], tempsCip[j + (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 * (int)1.000]);
		}
		for (int j = 0; j < 64; j++)
		{
			for (int k = 1; k < 1 * (int)1.000; k++)
			{
				evaluator.add_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j + (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 * k]);
			}
			for (int k = 1; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], k, galois_keys, tempCip);
				evaluator.add_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], tempCip);
			}
			evaluator.multiply_plain_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], one_filterPlain);
		}
		for (int j = 1; j < 64; j++)
		{
			evaluator.rotate_vector_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], -updatestep[j], galois_keys);
			evaluator.add_inplace(tempsCip[0], tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j]);
		}
		wb_4_deltaCip = tempsCip[0];
		evaluator.rescale_to_next_inplace(wb_4_deltaCip);
		//数据提供者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < hidenode; k++)
			{
				for (int m = 0; m < hidenode; m++)
				{
					w_4_delta[k][m] = 0;
				}
				b_4_delta[k] = 0;
			}
			for (int k = 0; k < slotCount; k++)
			{
				noises[j - 1][0][k] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
			}
			ckks_encoder.encode(noises[j - 1][0], levelFour, levelFourScale, tempPlain);
			decryptor[j]->partial_decrypt(wb_4_deltaCip, tempPlain, pDecPlain[0][j - 1]);
			for (int k = 0; k < hidenode; k++)
			{
				for (int m = 0; m < hidenode; m++)
				{
					w_4_delta[m][k] -= noises[j - 1][0][(long long)m * 256 % 4096 + updatestep[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 2 + m / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(m - k, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
				}
				b_4_delta[k] -= noises[j - 1][0][(long long)hidenode * 256 % 4096 + updatestep[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 2 + hidenode / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(hidenode - k, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
			}
			noisesCip[j - 1][0] = zeroCip;
			parameterEncode(noisesCip[j - 1][0], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w_4_delta, b_4_delta);
			tempCip = noisesCip[j - 1][0];
			evaluator.add_inplace(noisesCip[j - 1][0], wb_4_NoiseCip[j - 1]);
			evaluator.add_inplace(wb_4_NoiseCip[j - 1], tempCip);
			noisesCip[j - 1][1] = zeroCip;
			parameterEncodeTranspose(noisesCip[j - 1][1], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w_4_delta);
			tempCip = noisesCip[j - 1][1];
			evaluator.add_inplace(noisesCip[j - 1][1], wb_4_TranNoiseCip[j - 1]);
			evaluator.add_inplace(wb_4_TranNoiseCip[j - 1], tempCip);
		}
		//数据提供者：将2×(partyNum - 1)个噪声密文发送，partyNum - 1个部分解密发送
		//研究者：完全解密并重新编码
		decryptor[0]->final_decrypt(wb_4_deltaCip, pDecPlain[0], partyNum - 1, tempPlain);
		ckks_encoder.decode(tempPlain, noises[1][0]);
		for (int j = 0; j < hidenode; j++)
		{
			for (int k = 0; k < hidenode; k++)
			{
				w_4[k][j] += noises[1][0][(long long)k * 256 % 4096 + updatestep[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 2 + k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(k - j, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
			}
			b_4[j] += noises[1][0][(long long)hidenode * 256 % 4096 + updatestep[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 2 + hidenode / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(hidenode - j, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
		}
		wb_4Cip = zeroCip;
		parameterEncode(wb_4Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w_4, b_4);
		for (int j = 1; j < partyNum; j++)
		{
			evaluator.add_inplace(wb_4Cip, noisesCip[j - 1][0]);
		}
		wb_4TranCip = zeroCip;
		parameterEncodeTranspose(wb_4TranCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w_4);
		for (int j = 1; j < partyNum; j++)
		{
			evaluator.add_inplace(wb_4TranCip, noisesCip[j - 1][1]);
		}
		//数据提供者：保存隐含层之间的梯度(16倍数_3)
		Pos = 0;
		for (int j = 0; j < batchsize * time / 16; j++)
		{
			for (int k = 0; k < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k++)
			{
				coeff[Pos + j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) % 2][j % (int)((slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) - pow(2, ceil(log2(hidenode * hidenode / 4) / 2))) * 2) / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2 + 1)) * pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j % (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + (long long)k * slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2))] = (outcome[(i * batchsize + j / (time / 16)) % row] * 0.84 + 0.16) * alpha;
			}
			if (j % (int)((slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) - pow(2, ceil(log2(hidenode * hidenode / 4) / 2))) * 2) == (slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) - pow(2, ceil(log2(hidenode * hidenode / 4) / 2))) * 2 - 1)
			{
				Pos += 2;
			}
		}
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			ckks_encoder.encode(coeff[j], levelOne, levelOneScale, tempsPlain[j]);
		}
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 8 * (int)1.000; k++)
			{
				for (int m = 1; m < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); m++)
				{
					evaluator.rotate_vector(layer_1_3_deltaOriCip[j][k][m - 1], -slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), galois_keys, layer_1_3_deltaOriCip[j][k][m]);
				}
			}
		}
		for (int j = 0; j < 2 * (int)1.000; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				evaluator.multiply_plain(layer_1_3_0OriCip[k + j * 4], tempsPlain[k % 2 + j * 2], tempCip);
				evaluator.rescale_to_next_inplace(tempCip);
				for (int m = 0; m < 2; m++)
				{
					for (int n = 0; n < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); n++)
					{
						evaluator.multiply(layer_1_3_deltaOriCip[1][m * 2 + k % 2 + j * 4][n], tempCip, tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
						evaluator.relinearize_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n], relin_keys);
						evaluator.rescale_to_next_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
					}
				}
			}
		}
		for (int j = 0; j < 2 * (int)1.000; j++)
		{
			for (int k = 0; k < 64; k++)
			{
				evaluator.add_inplace(tempsCip[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8], tempsCip[(k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + 1) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8]);
			}
		}
		for (int j = 0; j < 64; j++)
		{
			for (int k = 1; k < 2 * (int)1.000; k++)
			{
				evaluator.add_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j + (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 * k]);
			}
			for (int k = 1; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], k, galois_keys, tempCip);
				evaluator.add_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], tempCip);
			}
			evaluator.multiply_plain_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], one_filterPlain);
		}
		for (int j = 1; j < 64; j++)
		{
			evaluator.rotate_vector_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], -updatestep[j], galois_keys);
			evaluator.add_inplace(tempsCip[0], tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j]);
		}
		wh_3_deltaCip = tempsCip[0];
		evaluator.rescale_to_next_inplace(wh_3_deltaCip);
		//数据提供者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < hidenode; k++)
			{
				for (int m = 0; m < hidenode; m++)
				{
					wh_3_delta[k][m] = 0;
				}
			}
			for (int k = 0; k < slotCount; k++)
			{
				noises[j - 1][0][k] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
			}
			ckks_encoder.encode(noises[j - 1][0], levelFour, levelFourScale, tempPlain);
			decryptor[j]->partial_decrypt(wh_3_deltaCip, tempPlain, pDecPlain[0][j - 1]);
			for (int k = 0; k < hidenode; k++)
			{
				for (int m = 0; m < hidenode; m++)
				{
					wh_3_delta[m][k] -= noises[j - 1][0][(long long)m * 256 % 4096 + updatestep[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 2 + m / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(m - k, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
				}
			}
			noisesCip[j - 1][0] = zeroCip;
			parameterEncode(noisesCip[j - 1][0], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_3_delta);
			tempCip = noisesCip[j - 1][0];
			evaluator.add_inplace(noisesCip[j - 1][0], wh_3_NoiseCip[j - 1]);
			evaluator.add_inplace(wh_3_NoiseCip[j - 1], tempCip);
			noisesCip[j - 1][1] = zeroCip;
			parameterEncodeTranspose(noisesCip[j - 1][1], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_3_delta);
			tempCip = noisesCip[j - 1][1];
			evaluator.add_inplace(noisesCip[j - 1][1], wh_3_TranNoiseCip[j - 1]);
			evaluator.add_inplace(wh_3_TranNoiseCip[j - 1], tempCip);
		}
		//数据提供者：将2×(partyNum - 1)个噪声密文发送，partyNum - 1个部分解密发送
		//研究者：完全解密并重新编码
		decryptor[0]->final_decrypt(wh_3_deltaCip, pDecPlain[0], partyNum - 1, tempPlain);
		ckks_encoder.decode(tempPlain, noises[1][0]);
		for (int j = 0; j < hidenode; j++)
		{
			for (int k = 0; k < hidenode; k++)
			{
				wh_3[k][j] += noises[1][0][(long long)k * 256 % 4096 + updatestep[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 2 + k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(k - j, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
			}
		}
		wh_3Cip = zeroCip;
		parameterEncode(wh_3Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_3);
		for (int j = 1; j < partyNum; j++)
		{
			evaluator.add_inplace(wh_3Cip, noisesCip[j - 1][0]);
		}
		wh_3TranCip = zeroCip;
		parameterEncodeTranspose(wh_3TranCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_3);
		for (int j = 1; j < partyNum; j++)
		{
			evaluator.add_inplace(wh_3TranCip, noisesCip[j - 1][1]);
		}
		//数据提供者：保存输入层和隐含层之间的梯度(16倍数_3)
		for (int j = 0; j < 2 * (int)1.000; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				evaluator.multiply_plain(layer_1_2_1OriCip[k + j * 4 + 8 * (int)1.000], tempsPlain[k % 2 + j * 2], tempCip);
				evaluator.rescale_to_next_inplace(tempCip);
				for (int m = 0; m < 2; m++)
				{
					for (int n = 0; n < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); n++)
					{
						evaluator.multiply(layer_1_3_deltaOriCip[1][m * 2 + k % 2 + j * 4][n], tempCip, tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
						evaluator.relinearize_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n], relin_keys);
						evaluator.rescale_to_next_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
					}
				}
			}
		}
		for (int j = 0; j < 2 * (int)1.000; j++)
		{
			for (int k = 0; k < 64; k++)
			{
				evaluator.add_inplace(tempsCip[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8], tempsCip[(k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + 1) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8]);
			}
		}
		for (int j = 0; j < (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 16 * (int)1.000; j++)
		{
			tempsCip[j + (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 16 * (int)1.000] = tempsCip[j];
		}
		//数据提供者：保存输入层和隐含层之间的梯度(16余8_3)
		for (int j = 0; j < 2 * (int)1.000; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				evaluator.multiply_plain(layer_1_2_1OriCip[k + j * 4], tempsPlain[k % 2 + j * 2], tempCip);
				evaluator.rescale_to_next_inplace(tempCip);
				for (int m = 0; m < 2; m++)
				{
					for (int n = 0; n < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); n++)
					{
						evaluator.multiply(layer_1_3_deltaOriCip[0][m * 2 + k % 2 + j * 4][n], tempCip, tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
						evaluator.relinearize_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n], relin_keys);
						evaluator.rescale_to_next_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
					}
				}
			}
		}
		for (int j = 0; j < 2 * (int)1.000; j++)
		{
			for (int k = 0; k < 64; k++)
			{
				evaluator.add_inplace(tempsCip[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8], tempsCip[(k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + 1) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8]);
			}
		}
		for (int j = 0; j < (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 16 * (int)1.000; j++)
		{
			evaluator.add_inplace(tempsCip[j], tempsCip[j + (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 16 * (int)1.000]);
		}
		for (int j = 0; j < 64; j++)
		{
			for (int k = 1; k < 2 * (int)1.000; k++)
			{
				evaluator.add_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j + (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 * k]);
			}
			for (int k = 1; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], k, galois_keys, tempCip);
				evaluator.add_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], tempCip);
			}
			evaluator.multiply_plain_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], one_filterPlain);
		}
		for (int j = 1; j < 64; j++)
		{
			evaluator.rotate_vector_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], -updatestep[j], galois_keys);
			evaluator.add_inplace(tempsCip[0], tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j]);
		}
		wb_3_deltaCip = tempsCip[0];
		evaluator.rescale_to_next_inplace(wb_3_deltaCip);
		//数据提供者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < hidenode; k++)
			{
				for (int m = 0; m < hidenode; m++)
				{
					w_3_delta[k][m] = 0;
				}
				b_3_delta[k] = 0;
			}
			for (int k = 0; k < slotCount; k++)
			{
				noises[j - 1][0][k] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
			}
			ckks_encoder.encode(noises[j - 1][0], levelFour, levelFourScale, tempPlain);
			decryptor[j]->partial_decrypt(wb_3_deltaCip, tempPlain, pDecPlain[0][j - 1]);
			for (int k = 0; k < hidenode; k++)
			{
				for (int m = 0; m < hidenode; m++)
				{
					w_3_delta[m][k] -= noises[j - 1][0][(long long)m * 256 % 4096 + updatestep[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 2 + m / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(m - k, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
				}
				b_3_delta[k] -= noises[j - 1][0][(long long)hidenode * 256 % 4096 + updatestep[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 2 + hidenode / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(hidenode - k, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
			}
			noisesCip[j - 1][0] = zeroCip;
			parameterEncode(noisesCip[j - 1][0], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w_3_delta, b_3_delta);
			tempCip = noisesCip[j - 1][0];
			evaluator.add_inplace(noisesCip[j - 1][0], wb_3_NoiseCip[j - 1]);
			evaluator.add_inplace(wb_3_NoiseCip[j - 1], tempCip);
			noisesCip[j - 1][1] = zeroCip;
			parameterEncodeTranspose(noisesCip[j - 1][1], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w_3_delta);
			tempCip = noisesCip[j - 1][1];
			evaluator.add_inplace(noisesCip[j - 1][1], wb_3_TranNoiseCip[j - 1]);
			evaluator.add_inplace(wb_3_TranNoiseCip[j - 1], tempCip);
		}
		//数据提供者：将2×(partyNum - 1)个噪声密文发送，partyNum - 1个部分解密发送
		//研究者：完全解密并重新编码
		decryptor[0]->final_decrypt(wb_3_deltaCip, pDecPlain[0], partyNum - 1, tempPlain);
		ckks_encoder.decode(tempPlain, noises[1][0]);
		for (int j = 0; j < hidenode; j++)
		{
			for (int k = 0; k < hidenode; k++)
			{
				w_3[k][j] += noises[1][0][(long long)k * 256 % 4096 + updatestep[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 2 + k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(k - j, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
			}
			b_3[j] += noises[1][0][(long long)hidenode * 256 % 4096 + updatestep[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 2 + hidenode / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(hidenode - j, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
		}
		wb_3Cip = zeroCip;
		parameterEncode(wb_3Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w_3, b_3);
		for (int j = 1; j < partyNum; j++)
		{
			evaluator.add_inplace(wb_3Cip, noisesCip[j - 1][0]);
		}
		wb_3TranCip = zeroCip;
		parameterEncodeTranspose(wb_3TranCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w_3);
		for (int j = 1; j < partyNum; j++)
		{
			evaluator.add_inplace(wb_3TranCip, noisesCip[j - 1][1]);
		}
		//数据提供者：保存隐含层之间的梯度(8倍数_2)
		Pos = 0;
		for (int j = 0; j < batchsize * time / 8; j++)
		{
			for (int k = 0; k < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k++)
			{
				coeff[Pos + j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) % 2][j % (int)((slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) - pow(2, ceil(log2(hidenode * hidenode / 4) / 2))) * 2) / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2 + 1)) * pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j % (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + (long long)k * slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2))] = (outcome[(i * batchsize + j / (time / 8)) % row] * 0.84 + 0.16) * alpha;
			}
			if (j % (int)((slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) - pow(2, ceil(log2(hidenode * hidenode / 4) / 2))) * 2) == (slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) - pow(2, ceil(log2(hidenode * hidenode / 4) / 2))) * 2 - 1)
			{
				Pos += 2;
			}
		}
		for (int j = 0; j < 8 * (int)1.000; j++)
		{
			ckks_encoder.encode(coeff[j], levelOne, levelOneScale, tempsPlain[j]);
		}
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 16 * (int)1.000; k++)
			{
				for (int m = 1; m < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); m++)
				{
					evaluator.rotate_vector(layer_1_2_deltaOriCip[j][k][m - 1], -slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), galois_keys, layer_1_2_deltaOriCip[j][k][m]);
				}
			}
		}
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				evaluator.multiply_plain(layer_1_2_0OriCip[k + j * 4], tempsPlain[k % 2 + j * 2], tempCip);
				evaluator.rescale_to_next_inplace(tempCip);
				for (int m = 0; m < 2; m++)
				{
					for (int n = 0; n < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); n++)
					{
						evaluator.multiply(layer_1_2_deltaOriCip[1][m * 2 + k % 2 + j * 4][n], tempCip, tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
						evaluator.relinearize_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n], relin_keys);
						evaluator.rescale_to_next_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
					}
				}
			}
		}
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			for (int k = 0; k < 64; k++)
			{
				evaluator.add_inplace(tempsCip[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8], tempsCip[(k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + 1) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8]);
			}
		}
		for (int j = 0; j < 64; j++)
		{
			for (int k = 1; k < 4 * (int)1.000; k++)
			{
				evaluator.add_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j + (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 * k]);
			}
			for (int k = 1; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], k, galois_keys, tempCip);
				evaluator.add_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], tempCip);
			}
			evaluator.multiply_plain_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], one_filterPlain);
		}
		for (int j = 1; j < 64; j++)
		{
			evaluator.rotate_vector_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], -updatestep[j], galois_keys);
			evaluator.add_inplace(tempsCip[0], tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j]);
		}
		wh_2_deltaCip = tempsCip[0];
		evaluator.rescale_to_next_inplace(wh_2_deltaCip);
		//数据提供者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < hidenode; k++)
			{
				for (int m = 0; m < hidenode; m++)
				{
					wh_2_delta[k][m] = 0;
				}
			}
			for (int k = 0; k < slotCount; k++)
			{
				noises[j - 1][0][k] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
			}
			ckks_encoder.encode(noises[j - 1][0], levelFour, levelFourScale, tempPlain);
			decryptor[j]->partial_decrypt(wh_2_deltaCip, tempPlain, pDecPlain[0][j - 1]);
			for (int k = 0; k < hidenode; k++)
			{
				for (int m = 0; m < hidenode; m++)
				{
					wh_2_delta[m][k] -= noises[j - 1][0][(long long)m * 256 % 4096 + updatestep[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 2 + m / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(m - k, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
				}
			}
			noisesCip[j - 1][0] = zeroCip;
			parameterEncode(noisesCip[j - 1][0], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_2_delta);
			tempCip = noisesCip[j - 1][0];
			evaluator.add_inplace(noisesCip[j - 1][0], wh_2_NoiseCip[j - 1]);
			evaluator.add_inplace(wh_2_NoiseCip[j - 1], tempCip);
			noisesCip[j - 1][1] = zeroCip;
			parameterEncodeTranspose(noisesCip[j - 1][1], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_2_delta);
			tempCip = noisesCip[j - 1][1];
			evaluator.add_inplace(noisesCip[j - 1][1], wh_2_TranNoiseCip[j - 1]);
			evaluator.add_inplace(wh_2_TranNoiseCip[j - 1], tempCip);
		}
		//数据提供者：将2×(partyNum - 1)个噪声密文发送，partyNum - 1个部分解密发送
		//研究者：完全解密并重新编码
		decryptor[0]->final_decrypt(wh_2_deltaCip, pDecPlain[0], partyNum - 1, tempPlain);
		ckks_encoder.decode(tempPlain, noises[1][0]);
		for (int j = 0; j < hidenode; j++)
		{
			for (int k = 0; k < hidenode; k++)
			{
				wh_2[k][j] += noises[1][0][(long long)k * 256 % 4096 + updatestep[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 2 + k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(k - j, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
			}
		}
		wh_2Cip = zeroCip;
		parameterEncode(wh_2Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_2);
		for (int j = 1; j < partyNum; j++)
		{
			evaluator.add_inplace(wh_2Cip, noisesCip[j - 1][0]);
		}
		wh_2TranCip = zeroCip;
		parameterEncodeTranspose(wh_2TranCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_2);
		for (int j = 1; j < partyNum; j++)
		{
			evaluator.add_inplace(wh_2TranCip, noisesCip[j - 1][1]);
		}
		//数据提供者：保存输入层和隐含层之间的梯度(8倍数_2)
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				evaluator.multiply_plain(layer_1_1_1OriCip[k + j * 4 + 16 * (int)1.000], tempsPlain[k % 2 + j * 2], tempCip);
				evaluator.rescale_to_next_inplace(tempCip);
				for (int m = 0; m < 2; m++)
				{
					for (int n = 0; n < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); n++)
					{
						evaluator.multiply(layer_1_2_deltaOriCip[1][m * 2 + k % 2 + j * 4][n], tempCip, tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
						evaluator.relinearize_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n], relin_keys);
						evaluator.rescale_to_next_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
					}
				}
			}
		}
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			for (int k = 0; k < 64; k++)
			{
				evaluator.add_inplace(tempsCip[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8], tempsCip[(k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + 1) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8]);
			}
		}
		for (int j = 0; j < (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 32 * (int)1.000; j++)
		{
			tempsCip[j + (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 32 * (int)1.000] = tempsCip[j];
		}
		//数据提供者：保存输入层和隐含层之间的梯度(8余4_2)
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				evaluator.multiply_plain(layer_1_1_1OriCip[k + j * 4], tempsPlain[k % 2 + j * 2], tempCip);
				evaluator.rescale_to_next_inplace(tempCip);
				for (int m = 0; m < 2; m++)
				{
					for (int n = 0; n < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); n++)
					{
						evaluator.multiply(layer_1_2_deltaOriCip[0][m * 2 + k % 2 + j * 4][n], tempCip, tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
						evaluator.relinearize_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n], relin_keys);
						evaluator.rescale_to_next_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
					}
				}
			}
		}
		for (int j = 0; j < 4 * (int)1.000; j++)
		{
			for (int k = 0; k < 64; k++)
			{
				evaluator.add_inplace(tempsCip[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8], tempsCip[(k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + 1) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8]);
			}
		}
		for (int j = 0; j < (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 32 * (int)1.000; j++)
		{
			evaluator.add_inplace(tempsCip[j], tempsCip[j + (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 32 * (int)1.000]);
		}
		for (int j = 0; j < 64; j++)
		{
			for (int k = 1; k < 4 * (int)1.000; k++)
			{
				evaluator.add_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j + (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 * k]);
			}
			for (int k = 1; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], k, galois_keys, tempCip);
				evaluator.add_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], tempCip);
			}
			evaluator.multiply_plain_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], one_filterPlain);
		}
		for (int j = 1; j < 64; j++)
		{
			evaluator.rotate_vector_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], -updatestep[j], galois_keys);
			evaluator.add_inplace(tempsCip[0], tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j]);
		}
		wb_2_deltaCip = tempsCip[0];
		evaluator.rescale_to_next_inplace(wb_2_deltaCip);
		//数据提供者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < hidenode; k++)
			{
				for (int m = 0; m < hidenode; m++)
				{
					w_2_delta[k][m] = 0;
				}
				b_2_delta[k] = 0;
			}
			for (int k = 0; k < slotCount; k++)
			{
				noises[j - 1][0][k] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
			}
			ckks_encoder.encode(noises[j - 1][0], levelFour, levelFourScale, tempPlain);
			decryptor[j]->partial_decrypt(wb_2_deltaCip, tempPlain, pDecPlain[0][j - 1]);
			for (int k = 0; k < hidenode; k++)
			{
				for (int m = 0; m < hidenode; m++)
				{
					w_2_delta[m][k] -= noises[j - 1][0][(long long)m * 256 % 4096 + updatestep[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 2 + m / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(m - k, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
				}
				b_2_delta[k] -= noises[j - 1][0][(long long)hidenode * 256 % 4096 + updatestep[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 2 + hidenode / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(hidenode - k, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
			}
			noisesCip[j - 1][0] = zeroCip;
			parameterEncode(noisesCip[j - 1][0], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w_2_delta, b_2_delta);
			tempCip = noisesCip[j - 1][0];
			evaluator.add_inplace(noisesCip[j - 1][0], wb_2_NoiseCip[j - 1]);
			evaluator.add_inplace(wb_2_NoiseCip[j - 1], tempCip);
			noisesCip[j - 1][1] = zeroCip;
			parameterEncodeTranspose(noisesCip[j - 1][1], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w_2_delta);
			tempCip = noisesCip[j - 1][1];
			evaluator.add_inplace(noisesCip[j - 1][1], wb_2_TranNoiseCip[j - 1]);
			evaluator.add_inplace(wb_2_TranNoiseCip[j - 1], tempCip);
		}
		//数据提供者：将2×(partyNum - 1)个噪声密文发送，partyNum - 1个部分解密发送
		//研究者：完全解密并重新编码
		decryptor[0]->final_decrypt(wb_2_deltaCip, pDecPlain[0], partyNum - 1, tempPlain);
		ckks_encoder.decode(tempPlain, noises[1][0]);
		for (int j = 0; j < hidenode; j++)
		{
			for (int k = 0; k < hidenode; k++)
			{
				w_2[k][j] += noises[1][0][(long long)k * 256 % 4096 + updatestep[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 2 + k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(k - j, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
			}
			b_2[j] += noises[1][0][(long long)hidenode * 256 % 4096 + updatestep[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 2 + hidenode / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(hidenode - j, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
		}
		wb_2Cip = zeroCip;
		parameterEncode(wb_2Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w_2, b_2);
		for (int j = 1; j < partyNum; j++)
		{
			evaluator.add_inplace(wb_2Cip, noisesCip[j - 1][0]);
		}
		wb_2TranCip = zeroCip;
		parameterEncodeTranspose(wb_2TranCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w_2);
		for (int j = 1; j < partyNum; j++)
		{
			evaluator.add_inplace(wb_2TranCip, noisesCip[j - 1][1]);
		}
		//数据提供者：保存隐含层之间的梯度(4倍数_1)
		Pos = 0;
		for (int j = 0; j < batchsize * time / 4; j++)
		{
			for (int k = 0; k < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k++)
			{
				coeff[Pos + j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) % 2][j % (int)((slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) - pow(2, ceil(log2(hidenode * hidenode / 4) / 2))) * 2) / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2 + 1)) * pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j % (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + (long long)k * slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2))] = (outcome[(i * batchsize + j / (time / 4)) % row] * 0.84 + 0.16) * alpha;
			}
			if (j % (int)((slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) - pow(2, ceil(log2(hidenode * hidenode / 4) / 2))) * 2) == (slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) - pow(2, ceil(log2(hidenode * hidenode / 4) / 2))) * 2 - 1)
			{
				Pos += 2;
			}
		}
		for (int j = 0; j < 14 * (int)1.000; j++)
		{
			ckks_encoder.encode(coeff[j], levelOne, levelOneScale, tempsPlain[j]);
		}
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 28 * (int)1.000; k++)
			{
				for (int m = 1; m < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); m++)
				{
					evaluator.rotate_vector(layer_1_1_deltaOriCip[j][k][m - 1], -slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), galois_keys, layer_1_1_deltaOriCip[j][k][m]);
				}
			}
		}
		for (int j = 0; j < 7 * (int)1.000; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				evaluator.multiply_plain(layer_1_1_0OriCip[k + j * 4], tempsPlain[k % 2 + j * 2], tempCip);
				evaluator.rescale_to_next_inplace(tempCip);
				for (int m = 0; m < 2; m++)
				{
					for (int n = 0; n < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); n++)
					{
						evaluator.multiply(layer_1_1_deltaOriCip[1][m * 2 + k % 2 + j * 4][n], tempCip, tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
						evaluator.relinearize_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n], relin_keys);
						evaluator.rescale_to_next_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
					}
				}
			}
		}
		for (int j = 0; j < 7 * (int)1.000; j++)
		{
			for (int k = 0; k < 64; k++)
			{
				evaluator.add_inplace(tempsCip[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8], tempsCip[(k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + 1) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8]);
			}
		}
		for (int j = 0; j < 64; j++)
		{
			for (int k = 1; k < 7 * (int)1.000; k++)
			{
				evaluator.add_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j + (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 * k]);
			}
			for (int k = 1; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], k, galois_keys, tempCip);
				evaluator.add_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], tempCip);
			}
			evaluator.multiply_plain_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], one_filterPlain);
		}
		for (int j = 1; j < 64; j++)
		{
			evaluator.rotate_vector_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], -updatestep[j], galois_keys);
			evaluator.add_inplace(tempsCip[0], tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j]);
		}
		wh_1_deltaCip = tempsCip[0];
		evaluator.rescale_to_next_inplace(wh_1_deltaCip);
		//数据提供者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < hidenode; k++)
			{
				for (int m = 0; m < hidenode; m++)
				{
					wh_1_delta[k][m] = 0;
				}
			}
			for (int k = 0; k < slotCount; k++)
			{
				noises[j - 1][0][k] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
			}
			ckks_encoder.encode(noises[j - 1][0], levelFour, levelFourScale, tempPlain);
			decryptor[j]->partial_decrypt(wh_1_deltaCip, tempPlain, pDecPlain[0][j - 1]);
			for (int k = 0; k < hidenode; k++)
			{
				for (int m = 0; m < hidenode; m++)
				{
					wh_1_delta[m][k] -= noises[j - 1][0][(long long)m * 256 % 4096 + updatestep[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 2 + m / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(m - k, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
				}
			}
			noisesCip[j - 1][0] = zeroCip;
			parameterEncode(noisesCip[j - 1][0], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_1_delta);
			tempCip = noisesCip[j - 1][0];
			evaluator.add_inplace(noisesCip[j - 1][0], wh_1_NoiseCip[j - 1]);
			evaluator.add_inplace(wh_1_NoiseCip[j - 1], tempCip);
			noisesCip[j - 1][1] = zeroCip;
			parameterEncodeTranspose(noisesCip[j - 1][1], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_1_delta);
			tempCip = noisesCip[j - 1][1];
			evaluator.add_inplace(noisesCip[j - 1][1], wh_1_TranNoiseCip[j - 1]);
			evaluator.add_inplace(wh_1_TranNoiseCip[j - 1], tempCip);
		}
		//数据提供者：将2×(partyNum - 1)个噪声密文发送，partyNum - 1个部分解密发送
		//研究者：完全解密并重新编码
		decryptor[0]->final_decrypt(wh_1_deltaCip, pDecPlain[0], partyNum - 1, tempPlain);
		ckks_encoder.decode(tempPlain, noises[1][0]);
		for (int j = 0; j < hidenode; j++)
		{
			for (int k = 0; k < hidenode; k++)
			{
				wh_1[k][j] += noises[1][0][(long long)k * 256 % 4096 + updatestep[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 2 + k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(k - j, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
			}
		}
		wh_1Cip = zeroCip;
		parameterEncode(wh_1Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_1);
		for (int j = 1; j < partyNum; j++)
		{
			evaluator.add_inplace(wh_1Cip, noisesCip[j - 1][0]);
		}
		wh_1TranCip = zeroCip;
		parameterEncodeTranspose(wh_1TranCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_1);
		for (int j = 1; j < partyNum; j++)
		{
			evaluator.add_inplace(wh_1TranCip, noisesCip[j - 1][1]);
		}
		//数据提供者：保存输入层和隐含层之间的梯度(4倍数_1)
		for (int j = 0; j < 7 * (int)1.000; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				evaluator.multiply_plain(layer_1_0_1OriCip[k + j * 4 + 28 * (int)1.000], tempsPlain[k % 2 + j * 2], tempCip);
				evaluator.rescale_to_next_inplace(tempCip);
				for (int m = 0; m < 2; m++)
				{
					for (int n = 0; n < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); n++)
					{
						evaluator.multiply(layer_1_1_deltaOriCip[1][m * 2 + k % 2 + j * 4][n], tempCip, tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
						evaluator.relinearize_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n], relin_keys);
						evaluator.rescale_to_next_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
					}
				}
			}
		}
		for (int j = 0; j < 7 * (int)1.000; j++)
		{
			for (int k = 0; k < 64; k++)
			{
				evaluator.add_inplace(tempsCip[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8], tempsCip[(k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + 1) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8]);
			}
		}
		for (int j = 0; j < (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 56 * (int)1.000; j++)
		{
			tempsCip[j + (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 56 * (int)1.000] = tempsCip[j];
		}
		//数据提供者：保存输入层和隐含层之间的梯度(4余2_1)
		for (int j = 0; j < 7 * (int)1.000; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				evaluator.multiply_plain(layer_1_0_1OriCip[k + j * 4], tempsPlain[k % 2 + j * 2], tempCip);
				evaluator.rescale_to_next_inplace(tempCip);
				for (int m = 0; m < 2; m++)
				{
					for (int n = 0; n < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); n++)
					{
						evaluator.multiply(layer_1_1_deltaOriCip[0][m * 2 + k % 2 + j * 4][n], tempCip, tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
						evaluator.relinearize_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n], relin_keys);
						evaluator.rescale_to_next_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
					}
				}
			}
		}
		for (int j = 0; j < 7 * (int)1.000; j++)
		{
			for (int k = 0; k < 64; k++)
			{
				evaluator.add_inplace(tempsCip[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8], tempsCip[(k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + 1) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8]);
			}
		}
		for (int j = 0; j < (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 56 * (int)1.000; j++)
		{
			evaluator.add_inplace(tempsCip[j], tempsCip[j + (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 56 * (int)1.000]);
		}
		for (int j = 0; j < 64; j++)
		{
			for (int k = 1; k < 7 * (int)1.000; k++)
			{
				evaluator.add_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j + (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 * k]);
			}
			for (int k = 1; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], k, galois_keys, tempCip);
				evaluator.add_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], tempCip);
			}
			evaluator.multiply_plain_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], one_filterPlain);
		}
		for (int j = 1; j < 64; j++)
		{
			evaluator.rotate_vector_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], -updatestep[j], galois_keys);
			evaluator.add_inplace(tempsCip[0], tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j]);
		}
		wb_1_deltaCip = tempsCip[0];
		evaluator.rescale_to_next_inplace(wb_1_deltaCip);
		//数据提供者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < hidenode; k++)
			{
				for (int m = 0; m < hidenode; m++)
				{
					w_1_delta[k][m] = 0;
				}
				b_1_delta[k] = 0;
			}
			for (int k = 0; k < slotCount; k++)
			{
				noises[j - 1][0][k] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
			}
			ckks_encoder.encode(noises[j - 1][0], levelFour, levelFourScale, tempPlain);
			decryptor[j]->partial_decrypt(wb_1_deltaCip, tempPlain, pDecPlain[0][j - 1]);
			for (int k = 0; k < hidenode; k++)
			{
				for (int m = 0; m < hidenode; m++)
				{
					w_1_delta[m][k] -= noises[j - 1][0][(long long)m * 256 % 4096 + updatestep[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 2 + m / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(m - k, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
				}
				b_1_delta[k] -= noises[j - 1][0][(long long)hidenode * 256 % 4096 + updatestep[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 2 + hidenode / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(hidenode - k, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
			}
			noisesCip[j - 1][0] = zeroCip;
			parameterEncode(noisesCip[j - 1][0], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w_1_delta, b_1_delta);
			tempCip = noisesCip[j - 1][0];
			evaluator.add_inplace(noisesCip[j - 1][0], wb_1_NoiseCip[j - 1]);
			evaluator.add_inplace(wb_1_NoiseCip[j - 1], tempCip);
			noisesCip[j - 1][1] = zeroCip;
			parameterEncodeTranspose(noisesCip[j - 1][1], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w_1_delta);
			tempCip = noisesCip[j - 1][1];
			evaluator.add_inplace(noisesCip[j - 1][1], wb_1_TranNoiseCip[j - 1]);
			evaluator.add_inplace(wb_1_TranNoiseCip[j - 1], tempCip);
		}
		//数据提供者：将2×(partyNum - 1)个噪声密文发送，partyNum - 1个部分解密发送
		//研究者：完全解密并重新编码
		decryptor[0]->final_decrypt(wb_1_deltaCip, pDecPlain[0], partyNum - 1, tempPlain);
		ckks_encoder.decode(tempPlain, noises[1][0]);
		for (int j = 0; j < hidenode; j++)
		{
			for (int k = 0; k < hidenode; k++)
			{
				w_1[k][j] += noises[1][0][(long long)k * 256 % 4096 + updatestep[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 2 + k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(k - j, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
			}
			b_1[j] += noises[1][0][(long long)hidenode * 256 % 4096 + updatestep[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 2 + hidenode / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(hidenode - j, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
		}
		wb_1Cip = zeroCip;
		parameterEncode(wb_1Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w_1, b_1);
		for (int j = 1; j < partyNum; j++)
		{
			evaluator.add_inplace(wb_1Cip, noisesCip[j - 1][0]);
		}
		wb_1TranCip = zeroCip;
		parameterEncodeTranspose(wb_1TranCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, w_1);
		for (int j = 1; j < partyNum; j++)
		{
			evaluator.add_inplace(wb_1TranCip, noisesCip[j - 1][1]);
		}
		//数据提供者：保存隐含层之间的梯度(2倍数_0)
		Pos = 0;
		for (int j = 0; j < batchsize * time / 2; j++)
		{
			for (int k = 0; k < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k++)
			{
				coeff[Pos + j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) % 2][j % (int)((slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) - pow(2, ceil(log2(hidenode * hidenode / 4) / 2))) * 2) / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2 + 1)) * pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j % (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + (long long)k * slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2))] = (outcome[(i * batchsize + j / (time / 2)) % row] * 0.84 + 0.16) * alpha;
			}
			if (j % (int)((slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) - pow(2, ceil(log2(hidenode * hidenode / 4) / 2))) * 2) == (slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) - pow(2, ceil(log2(hidenode * hidenode / 4) / 2))) * 2 - 1)
			{
				Pos += 2;
			}
		}
		for (int j = 0; j < 26 * (int)1.000; j++)
		{
			ckks_encoder.encode(coeff[j], levelOne, levelOneScale, tempsPlain[j]);
		}
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 52 * (int)1.000; k++)
			{
				for (int m = 1; m < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); m++)
				{
					evaluator.rotate_vector(layer_1_0_deltaOriCip[j][k][m - 1], -slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), galois_keys, layer_1_0_deltaOriCip[j][k][m]);
				}
			}
		}
		for (int j = 0; j < 13 * (int)1.000; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				evaluator.multiply_plain(layer_1_0_0OriCip[k + j * 4], tempsPlain[k % 2 + j * 2], tempCip);
				evaluator.rescale_to_next_inplace(tempCip);
				for (int m = 0; m < 2; m++)
				{
					for (int n = 0; n < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); n++)
					{
						evaluator.multiply(layer_1_0_deltaOriCip[1][m * 2 + k % 2 + j * 4][n], tempCip, tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
						evaluator.relinearize_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n], relin_keys);
						evaluator.rescale_to_next_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 4 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 + n]);
					}
				}
			}
		}
		for (int j = 0; j < 13 * (int)1.000; j++)
		{
			for (int k = 0; k < 64; k++)
			{
				evaluator.add_inplace(tempsCip[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8], tempsCip[(k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + 1) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8]);
			}
		}
		for (int j = 0; j < 64; j++)
		{
			for (int k = 1; k < 13 * (int)1.000; k++)
			{
				evaluator.add_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j + (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 8 * k]);
			}
			for (int k = 1; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], k, galois_keys, tempCip);
				evaluator.add_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], tempCip);
			}
			evaluator.multiply_plain_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], one_filterPlain);
		}
		for (int j = 1; j < 64; j++)
		{
			evaluator.rotate_vector_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], -updatestep[j], galois_keys);
			evaluator.add_inplace(tempsCip[0], tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j]);
		}
		wh_0_deltaCip = tempsCip[0];
		evaluator.rescale_to_next_inplace(wh_0_deltaCip);
		//数据提供者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < hidenode; k++)
			{
				for (int m = 0; m < hidenode; m++)
				{
					wh_0_delta[k][m] = 0;
				}
			}
			for (int k = 0; k < slotCount; k++)
			{
				noises[j - 1][0][k] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
			}
			ckks_encoder.encode(noises[j - 1][0], levelFour, levelFourScale, tempPlain);
			decryptor[j]->partial_decrypt(wh_0_deltaCip, tempPlain, pDecPlain[0][j - 1]);
			for (int k = 0; k < hidenode; k++)
			{
				for (int m = 0; m < hidenode; m++)
				{
					wh_0_delta[m][k] -= noises[j - 1][0][(long long)m * 256 % 4096 + updatestep[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 2 + m / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(m - k, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
				}
			}
			noisesCip[j - 1][0] = zeroCip;
			parameterEncode(noisesCip[j - 1][0], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_0_delta);
			tempCip = noisesCip[j - 1][0];
			evaluator.add_inplace(noisesCip[j - 1][0], wh_0_NoiseCip[j - 1]);
			evaluator.add_inplace(wh_0_NoiseCip[j - 1], tempCip);
			noisesCip[j - 1][1] = zeroCip;
			parameterEncodeTranspose(noisesCip[j - 1][1], hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_0_delta);
			tempCip = noisesCip[j - 1][1];
			evaluator.add_inplace(noisesCip[j - 1][1], wh_0_TranNoiseCip[j - 1]);
			evaluator.add_inplace(wh_0_TranNoiseCip[j - 1], tempCip);
		}
		//数据提供者：将2×(partyNum - 1)个噪声密文发送，partyNum - 1个部分解密发送
		//研究者：完全解密并重新编码
		decryptor[0]->final_decrypt(wh_0_deltaCip, pDecPlain[0], partyNum - 1, tempPlain);
		ckks_encoder.decode(tempPlain, noises[1][0]);
		for (int j = 0; j < hidenode; j++)
		{
			for (int k = 0; k < hidenode; k++)
			{
				wh_0[k][j] += noises[1][0][(long long)k * 256 % 4096 + updatestep[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 2 + k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(k - j, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
			}
		}
		wh_0Cip = zeroCip;
		parameterEncode(wh_0Cip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_0);
		for (int j = 1; j < partyNum; j++)
		{
			evaluator.add_inplace(wh_0Cip, noisesCip[j - 1][0]);
		}
		wh_0TranCip = zeroCip;
		parameterEncodeTranspose(wh_0TranCip, hidenode, hidenode, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)), evaluator, ckks_encoder, empty, wh_0);
		for (int j = 1; j < partyNum; j++)
		{
			evaluator.add_inplace(wh_0TranCip, noisesCip[j - 1][1]);
		}
		//数据提供者：保存输入层和隐含层之间的梯度(2倍数_0)
		dataCode(layer_0_0_1, hidenode, col, pow(2, ceil(log2(hidenode * col / 6) / 2)), data, 1, i);
		for (int j = 0; j < 78 * (int)1.000; j++)
		{
			ckks_encoder.encode(layer_0_0_1[j], levelOne, levelOneScale, tempPlain);
			encryptor.encrypt(tempPlain, layer_0_0Cip[j]);
		}
		for (int j = 0; j < 13 * (int)1.000; j++)
		{
			for (int k = 0; k < 6; k++)
			{
				evaluator.multiply_plain(layer_0_0Cip[k + j * 6], tempsPlain[k % 2 + j * 2], tempCip);
				evaluator.rescale_to_next_inplace(tempCip);
				for (int m = 0; m < 2; m++)
				{
					for (int n = 0; n < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); n++)
					{
						evaluator.multiply(layer_1_0_deltaOriCip[1][m * 2 + k % 2 + j * 4][n], tempCip, tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 6 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 12 + n]);
						evaluator.relinearize_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 6 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 12 + n], relin_keys);
						evaluator.rescale_to_next_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 6 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 12 + n]);
					}
				}
			}
		}
		for (int j = 0; j < 13 * (int)1.000; j++)
		{
			for (int k = 0; k < 96; k++)
			{
				evaluator.add_inplace(tempsCip[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 12], tempsCip[(k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + 1) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 12]);
			}
		}
		for (int j = 0; j < (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 156 * (int)1.000; j++)
		{
			tempsCip[j + (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 156 * (int)1.000] = tempsCip[j];
		}
		//数据提供者：保存输入层和隐含层之间的梯度(2余1_0)
		dataCode(layer_0_0_0, hidenode, col, pow(2, ceil(log2(hidenode * col / 6) / 2)), data, 0, i);
		for (int j = 0; j < 78 * (int)1.000; j++)
		{
			ckks_encoder.encode(layer_0_0_0[j], levelOne, levelOneScale, tempPlain);
			encryptor.encrypt(tempPlain, layer_0_0Cip[j]);
		}
		for (int j = 0; j < 13 * (int)1.000; j++)
		{
			for (int k = 0; k < 6; k++)
			{
				evaluator.multiply_plain(layer_0_0Cip[k + j * 6], tempsPlain[k % 2 + j * 2], tempCip);
				evaluator.rescale_to_next_inplace(tempCip);
				for (int m = 0; m < 2; m++)
				{
					for (int n = 0; n < pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); n++)
					{
						evaluator.multiply(layer_1_0_deltaOriCip[0][m * 2 + k % 2 + j * 4][n], tempCip, tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 6 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 12 + n]);
						evaluator.relinearize_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 6 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 12 + n], relin_keys);
						evaluator.rescale_to_next_inplace(tempsCip[k * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + m * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 6 + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 12 + n]);
					}
				}
			}
		}
		for (int j = 0; j < 13 * (int)1.000; j++)
		{
			for (int k = 0; k < 96; k++)
			{
				evaluator.add_inplace(tempsCip[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 12], tempsCip[(k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + 1) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + k + j * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 12]);
			}
		}
		for (int j = 0; j < (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 156 * (int)1.000; j++)
		{
			evaluator.add_inplace(tempsCip[j], tempsCip[j + (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 156 * (int)1.000]);
		}
		for (int j = 0; j < 96; j++)
		{
			for (int k = 1; k < 13 * (int)1.000; k++)
			{
				evaluator.add_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j + (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 12 * k]);
			}
			for (int k = 1; k < slotCount / pow(2, ceil(log2(hidenode * hidenode / 4) / 2)); k *= 2)
			{
				evaluator.rotate_vector(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], k, galois_keys, tempCip);
				evaluator.add_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], tempCip);
			}
			evaluator.multiply_plain_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], one_filterPlain);
		}
		for (int j = 1; j < 96; j++)
		{
			evaluator.rotate_vector_inplace(tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j], -updatezero[j], galois_keys);
			evaluator.add_inplace(tempsCip[0], tempsCip[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + j]);
		}
		wb_0_deltaCip = tempsCip[0];
		evaluator.rescale_to_next_inplace(wb_0_deltaCip);
		//数据提供者：增加噪声并部分解密
		for (int j = 1; j < partyNum; j++)
		{
			for (int k = 0; k < hidenode; k++)
			{
				for (int m = 0; m < col; m++)
				{
					w_0_delta[m][k] = 0;
				}
				b_0_delta[k] = 0;
			}
			for (int k = 0; k < slotCount; k++)
			{
				noises[j - 1][0][k] = ((rand() << 15) + rand()) % 40000000 / 100.0 - 200000;
			}
			ckks_encoder.encode(noises[j - 1][0], levelFour, levelFourScale, tempPlain);
			decryptor[j]->partial_decrypt(wb_0_deltaCip, tempPlain, pDecPlain[0][j - 1]);
			for (int k = 0; k < hidenode; k++)
			{
				for (int m = 0; m < col; m++)
				{
					w_0_delta[m][k] -= noises[j - 1][0][(long long)m * 256 % 4096 + updatezero[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 3 + m / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(m - k, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
				}
				b_0_delta[k] -= noises[j - 1][0][(long long)col * 256 % 4096 + updatezero[k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 3 + col / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(col - k, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
			}
			noisesCip[j - 1][0] = zeroCip;
			parameterEncode(noisesCip[j - 1][0], hidenode, col, pow(2, ceil(log2(hidenode * col / 6) / 2)), evaluator, ckks_encoder, empty, w_0_delta, b_0_delta);
			tempCip = noisesCip[j - 1][0];
			evaluator.add_inplace(noisesCip[j - 1][0], wb_0_NoiseCip[j - 1]);
			evaluator.add_inplace(wb_0_NoiseCip[j - 1], tempCip);
		}
		//数据提供者：将partyNum - 1个噪声密文发送，partyNum - 1个部分解密发送
		//研究者：完全解密并重新编码
		decryptor[0]->final_decrypt(wb_0_deltaCip, pDecPlain[0], partyNum - 1, tempPlain);
		ckks_encoder.decode(tempPlain, noises[1][0]);
		for (int j = 0; j < hidenode; j++)
		{
			for (int k = 0; k < col; k++)
			{
				w_0[k][j] += noises[1][0][(long long)k * 256 % 4096 + updatezero[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 3 + k / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(k - j, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
			}
			b_0[j] += noises[1][0][(long long)col * 256 % 4096 + updatezero[j / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * 3 + col / (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) * (int)pow(2, ceil(log2(hidenode * hidenode / 4) / 2)) + mod(col - j, pow(2, ceil(log2(hidenode * hidenode / 4) / 2)))]] / pow(2, 15.5);
		}
		wb_0Cip = zeroCip;
		parameterEncode(wb_0Cip, hidenode, col, pow(2, ceil(log2(hidenode * col / 6) / 2)), evaluator, ckks_encoder, empty, w_0, b_0);
		for (int j = 1; j < partyNum; j++)
		{
			evaluator.add_inplace(wb_0Cip, noisesCip[j - 1][0]);
		}
		cout << "婷婷" << endl;

		if (i == 32 * 17 - 1)
		{
			Decryptor ceshijiemi(context, keygen.secret_key());
			Plaintext ceshiplain;
			static vector<double> ceshi[11];
			ceshijiemi.decrypt(wb_0Cip, ceshiplain);
			ckks_encoder.decode(ceshiplain, ceshi[0]);
			ceshijiemi.decrypt(wb_1Cip, ceshiplain);
			ckks_encoder.decode(ceshiplain, ceshi[1]);
			ceshijiemi.decrypt(wb_2Cip, ceshiplain);
			ckks_encoder.decode(ceshiplain, ceshi[2]);
			ceshijiemi.decrypt(wb_3Cip, ceshiplain);
			ckks_encoder.decode(ceshiplain, ceshi[3]);
			ceshijiemi.decrypt(wb_4Cip, ceshiplain);
			ckks_encoder.decode(ceshiplain, ceshi[4]);
			ceshijiemi.decrypt(wh_0Cip, ceshiplain);
			ckks_encoder.decode(ceshiplain, ceshi[5]);
			ceshijiemi.decrypt(wh_1Cip, ceshiplain);
			ckks_encoder.decode(ceshiplain, ceshi[6]);
			ceshijiemi.decrypt(wh_2Cip, ceshiplain);
			ckks_encoder.decode(ceshiplain, ceshi[7]);
			ceshijiemi.decrypt(wh_3Cip, ceshiplain);
			ckks_encoder.decode(ceshiplain, ceshi[8]);
			ceshijiemi.decrypt(wh_4Cip, ceshiplain);
			ckks_encoder.decode(ceshiplain, ceshi[9]);
			ceshijiemi.decrypt(wb1_Cip, ceshiplain);
			ckks_encoder.decode(ceshiplain, ceshi[10]);
			for (int j = 0; j < 10; j++)
			{
				for (int k = 0; k < 4096; k++)
				{
					ceshi[j][k] /= 3.735111370627234931;
				}
			}
			ofstream outceshi("D:\\Ceshijieguo\\error_" + to_string(i + 1) + ".txt");
			for (int j = 0; j < 11; j++)
			{
				for (int k = 0; k < 4096; k++)
				{
					outceshi << fixed << setprecision(15) << ceshi[j][k] << endl;
				}
			}
			outceshi.close();
			cout << "芳芳" << endl;
		}
	}
}