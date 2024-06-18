#include "pch.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <malloc.h>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

//Settings  
#define row 4000
#define col 36
#define timestep 48
#define hidenode 31

double activation(double x)
{
	return 0.0042 * x * x + 0.5 * x;
	//return 0.5 * (x + abs(x));
}

double dactivate(double x)
{
	return 0.0042 * 2 * x + 0.5;
	/*if (x < 0)
		return 0;
	else
		return 1;*/
}

double outactive(double x)
{
	if ((x) < -10)
		return 4.5397868702434395e-05;
	else if ((x) < -9)
		return 7.7996707283797329e-05 * x + 0.00082536494154040769;
	else if ((x) < -8)
		return 0.00021195555448024638 * x + 0.0020309945663084493;
	else if ((x) < -7)
		return 0.00057570106393416728 * x + 0.0049409586419398160;
	else if ((x) < -6)
		return 0.0015615719622341291 * x + 0.011842054930039550;
	else if ((x) < -5)
		return 0.0042202277676500807 * x + 0.027793989762535259;
	else if ((x) < -4)
		return 0.011293359037806703 * x + 0.063159646113318368;
	else if ((x) < -3)
		return 0.029439663215475222 * x + 0.13574486282399245;
	else if ((x) < -2)
		return 0.071777048844550773 * x + 0.26275701971121912;
	else if ((x) < -1)
		return 0.14973849934787756 * x + 0.41867992071787263;
	else if ((x) < 0)
		return 0.23105857863000490 * x + 0.5;
	else if ((x) < 1)
		return 0.23105857863000490 * x + 0.5;
	else if ((x) < 2)
		return 0.14973849934787742 * x + 0.58132007928212748;
	else if ((x) < 3)
		return 0.071777048844551050 * x + 0.73724298028878021;
	else if ((x) < 4)
		return 0.029439663215475087 * x + 0.86425513717600810;
	else if ((x) < 5)
		return 0.011293359037806816 * x + 0.93684035388668119;
	else if ((x) < 6)
		return 0.0042202277676500755 * x + 0.97220601023746489;
	else if ((x) < 7)
		return 0.0015615719622340540 * x + 0.98815794506996102;
	else if ((x) < 8)
		return 0.00057570106393423082 * x + 0.99505904135805978;
	else if ((x) < 9)
		return 0.00021195555448005887 * x + 0.99796900543369316;
	else if ((x) < 10)
		return 7.7996707283922717e-05 * x + 0.99917463505845838;
	else
		return 0.99995460213129761;
}

double doutactive(double x)
{
	return x * (1 - x);
}

double gaussrand()
{
	static double V1, V2, S;
	static int phase = 0;
	double X;
	if (phase == 0) {
		do {
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;
			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
		} while (S >= 1 || S == 0);
		X = V1 * sqrt(-2 * log(S) / S);
	}
	else
	{
		X = V2 * sqrt(-2 * log(S) / S);
	}
	phase = 1 - phase;
	return X;
}

void test(double average[col], double stdev[col], double b_0[hidenode], double b_1[hidenode], double b_2[hidenode], double b_3[hidenode], double b_4[hidenode], double w_0[col][hidenode],
	double w_1[hidenode][hidenode], double w_2[hidenode][hidenode], double w_3[hidenode][hidenode], double w_4[hidenode][hidenode], double wh_0[hidenode][hidenode], double wh_1[hidenode][hidenode],
	double wh_2[hidenode][hidenode], double wh_3[hidenode][hidenode], double wh_4[hidenode][hidenode], double b1, double w1[hidenode], double(&result)[4000][1000], int index, int loop, double& e)
{
	ifstream in3("D:\\TestTestData" + to_string(loop) + ".csv");
	static double testData[4000][col][timestep];
	int d1 = 0, a1 = 0;
	string line;
	getline(in3, line);
	while (getline(in3, line))
	{
		for (int i = 0; i < timestep; i++)
		{
			for (int j = 0; j < col; j++)
			{
				a1 = line.find(',');
				if (line.substr(0, a1) == "男_1" || line.substr(0, a1) == "男")
				{
					testData[d1][j][i] = 1;
				}
				else if (line.substr(0, a1) == "女_2" || line.substr(0, a1) == "女")
				{
					testData[d1][j][i] = 0;
				}
				else if (line.substr(0, a1) != "")
				{
					testData[d1][j][i] = atof(line.substr(0, a1).c_str());
				}
				else
				{
					testData[d1][j][i] = 1000000;
				}
				line.erase(0, (long long)a1 + 1);
			}
		}
		d1++;
	}
	in3.close();
	for (int i = 0; i < 4000; i++)
	{
		for (int j = 0; j < col; j++)
		{
			for (int k = 0; k < timestep; k++)
			{
				if (testData[i][j][k] != 1000000)
				{
					testData[i][j][k] = (testData[i][j][k] - average[j]) / stdev[j];
				}
				else
				{
					testData[i][j][k] = 0;
				}
			}
		}
	}
	ifstream in4("D:\\OutOutcome" + to_string(loop) + ".csv");
	static double outcome[4000];
	getline(in4, line);
	for (int i = 0; i < 4000; i++)
	{
		getline(in4, line);
		outcome[i] = atof(line.c_str());
	}
	in4.close();
	static double layer_0[col];
	static double layer_1_0[timestep][hidenode];
	static double layer_1_1[timestep / 2][hidenode];
	static double layer_1_2[timestep / 4][hidenode];
	static double layer_1_3[timestep / 8][hidenode];
	static double layer_1_4[timestep / 16][hidenode];
	double layer_2;
	//loss
	e = 0;
	for (int i = 0; i < 4000; i++)
	{
		//1,3,5,7,...,47
		for (int j = 0; j < timestep / 2; j++)
		{
			for (int k = 0; k < col; k++)
			{
				layer_0[k] = testData[i][k][2 * j];
			}
			for (int k = 0; k < hidenode; k++)
			{
				double o1 = b_0[k];
				for (int m = 0; m < col; m++)
				{
					o1 += layer_0[m] * w_0[m][k];
				}
				layer_1_0[2 * j][k] = activation(o1);
			}
		}
		//2,4,6,8,...,48
		for (int j = 0; j < timestep / 2; j++)
		{
			for (int k = 0; k < col; k++)
			{
				layer_0[k] = testData[i][k][2 * j + 1];
			}
			for (int k = 0; k < hidenode; k++)
			{
				double o1 = b_0[k];
				for (int m = 0; m < col; m++)
				{
					o1 += layer_0[m] * w_0[m][k];
				}
				for (int m = 0; m < hidenode; m++)
				{
					o1 += layer_1_0[2 * j][m] * wh_0[m][k];
				}
				layer_1_0[2 * j + 1][k] = activation(o1);
			}
		}
		//2,6,10,14,...,46
		for (int j = 0; j < timestep / 4; j++)
		{
			for (int k = 0; k < hidenode; k++)
			{
				double o1 = b_1[k];
				for (int m = 0; m < hidenode; m++)
				{
					o1 += layer_1_0[4 * j + 1][m] * w_1[m][k];
				}
				layer_1_1[2 * j][k] = activation(o1);
			}
		}
		//4,8,12,16,...,48
		for (int j = 0; j < timestep / 4; j++)
		{
			for (int k = 0; k < hidenode; k++)
			{
				double o1 = b_1[k];
				for (int m = 0; m < hidenode; m++)
				{
					o1 += layer_1_0[4 * j + 3][m] * w_1[m][k];
				}
				for (int m = 0; m < hidenode; m++)
				{
					o1 += layer_1_1[2 * j][m] * wh_1[m][k];
				}
				layer_1_1[2 * j + 1][k] = activation(o1);
			}
		}
		//4,12,20,28,36,44
		for (int j = 0; j < timestep / 8; j++)
		{
			for (int k = 0; k < hidenode; k++)
			{
				double o1 = b_2[k];
				for (int m = 0; m < hidenode; m++)
				{
					o1 += layer_1_1[4 * j + 1][m] * w_2[m][k];
				}
				layer_1_2[2 * j][k] = activation(o1);
			}
		}
		//8,16,24,32,40,48
		for (int j = 0; j < timestep / 8; j++)
		{
			for (int k = 0; k < hidenode; k++)
			{
				double o1 = b_2[k];
				for (int m = 0; m < hidenode; m++)
				{
					o1 += layer_1_1[4 * j + 3][m] * w_2[m][k];
				}
				for (int m = 0; m < hidenode; m++)
				{
					o1 += layer_1_2[2 * j][m] * wh_2[m][k];
				}
				layer_1_2[2 * j + 1][k] = activation(o1);
			}
		}
		//8,24,40
		for (int j = 0; j < timestep / 16; j++)
		{
			for (int k = 0; k < hidenode; k++)
			{
				double o1 = b_3[k];
				for (int m = 0; m < hidenode; m++)
				{
					o1 += layer_1_2[4 * j + 1][m] * w_3[m][k];
				}
				layer_1_3[2 * j][k] = activation(o1);
			}
		}
		//16,32,48
		for (int j = 0; j < timestep / 16; j++)
		{
			for (int k = 0; k < hidenode; k++)
			{
				double o1 = b_3[k];
				for (int m = 0; m < hidenode; m++)
				{
					o1 += layer_1_2[4 * j + 3][m] * w_3[m][k];
				}
				for (int m = 0; m < hidenode; m++)
				{
					o1 += layer_1_3[2 * j][m] * wh_3[m][k];
				}
				layer_1_3[2 * j + 1][k] = activation(o1);
			}
		}
		//16
		for (int j = 0; j < hidenode; j++)
		{
			double o1 = b_4[j];
			for (int k = 0; k < hidenode; k++)
			{
				o1 += layer_1_3[1][k] * w_4[k][j];
			}
			layer_1_4[0][j] = activation(o1);
		}
		//32
		for (int j = 0; j < hidenode; j++)
		{
			double o1 = b_4[j];
			for (int k = 0; k < hidenode; k++)
			{
				o1 += layer_1_3[3][k] * w_4[k][j];
			}
			for (int k = 0; k < hidenode; k++)
			{
				o1 += layer_1_4[0][k] * wh_4[k][j];
			}
			layer_1_4[1][j] = activation(o1);
		}
		//48
		for (int j = 0; j < hidenode; j++)
		{
			double o1 = b_4[j];
			for (int k = 0; k < hidenode; k++)
			{
				o1 += layer_1_3[5][k] * w_4[k][j];
			}
			for (int k = 0; k < hidenode; k++)
			{
				o1 += layer_1_4[1][k] * wh_4[k][j];
			}
			layer_1_4[2][j] = activation(o1);
		}
		//Output layer
		double o2 = b1;
		for (int j = 0; j < hidenode; j++)
		{
			o2 += layer_1_4[2][j] * w1[j];
		}
		layer_2 = outactive(o2);
		result[i][index] = layer_2;
		//loss computing
		e += -(outcome[i] * log(layer_2) + (1 - outcome[i]) * log(1 - layer_2)) * (0.145475 + 0.854525 * outcome[i]);
	}
}

int main(int argc, char* argv[])
{
	for (int loop = 0; loop < atoi(argv[1]) + 1; loop++)
	{
		ofstream outout("D:\\Convergence_SRNN" + to_string(loop) + "_" + to_string(atof(argv[2])) + ".csv");
		//data preparation
		ifstream in1("D:\\TestData" + to_string(loop) + ".csv");
		static double data[row][col][timestep];
		string line;
		size_t a1, d1 = 0;
		getline(in1, line);
		while (getline(in1, line))
		{
			for (int i = 0; i < timestep; i++)
			{
				for (int j = 0; j < col; j++)
				{
					a1 = line.find(',');
					if (line.substr(0, a1) == "男_1" || line.substr(0, a1) == "男")
					{
						data[d1][j][i] = 1;
					}
					else if (line.substr(0, a1) == "女_2" || line.substr(0, a1) == "女")
					{
						data[d1][j][i] = 0;
					}
					else if (line.substr(0, a1) != "")
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
			average[i] = 0;
			avect[i] = 0;
		}
		for (int i = 0; i < col; i++)
		{
			for (int j = 0; j < row; j++)
			{
				for (int k = 0; k < timestep; k++)
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
			stdev[i] = 0;
		}
		for (int i = 0; i < col; i++)
		{
			for (int j = 0; j < row; j++)
			{
				for (int k = 0; k < timestep; k++)
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
				for (int k = 0; k < timestep; k++)
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
		ifstream in2("D:\\Outcome" + to_string(loop) + ".csv");
		static double outcome[row];
		getline(in2, line);
		for (int i = 0; i < row; i++)
		{
			getline(in2, line);
			outcome[i] = atof(line.c_str());
		}
		in2.close();
		static double result[4000][1000];
		//training preparation
		int iter = row * 200;
		int minibatchsize = 250;
		static double b_0[hidenode], b_1[hidenode], b_2[hidenode], b_3[hidenode], b_4[hidenode];
		static double w_0[col][hidenode], w_1[hidenode][hidenode], w_2[hidenode][hidenode], w_3[hidenode][hidenode], w_4[hidenode][hidenode];
		static double wh_0[hidenode][hidenode], wh_1[hidenode][hidenode], wh_2[hidenode][hidenode], wh_3[hidenode][hidenode], wh_4[hidenode][hidenode];
		static double b1;
		static double w1[hidenode];
		static double layer_0[col];
		static double layer_1_0[timestep][hidenode];
		static double layer_1_1[timestep / 2][hidenode];
		static double layer_1_2[timestep / 4][hidenode];
		static double layer_1_3[timestep / 8][hidenode];
		static double layer_1_4[timestep / 16][hidenode];
		double layer_2;
		double layer_2_delta;
		static double layer_1_4_delta[timestep / 16][hidenode];
		static double layer_1_3_delta[timestep / 8][hidenode];
		static double layer_1_2_delta[timestep / 4][hidenode];
		static double layer_1_1_delta[timestep / 2][hidenode];
		static double layer_1_0_delta[timestep][hidenode];
		double b1_delta;
		static double w1_delta[hidenode];
		static double wh_0_delta[hidenode][hidenode], wh_1_delta[hidenode][hidenode], wh_2_delta[hidenode][hidenode], wh_3_delta[hidenode][hidenode], wh_4_delta[hidenode][hidenode];
		static double b_0_delta[hidenode], b_1_delta[hidenode], b_2_delta[hidenode], b_3_delta[hidenode], b_4_delta[hidenode];
		static double w_0_delta[col][hidenode], w_1_delta[hidenode][hidenode], w_2_delta[hidenode][hidenode], w_3_delta[hidenode][hidenode], w_4_delta[hidenode][hidenode];
		double alpha;
		static double o1_0[timestep][hidenode];
		static double o1_1[timestep / 2][hidenode];
		static double o1_2[timestep / 4][hidenode];
		static double o1_3[timestep / 8][hidenode];
		static double o1_4[timestep / 16][hidenode];
		double b1_delta_all = 0;
		static double w1_delta_all[hidenode];
		static double wh_0_delta_all[hidenode][hidenode], wh_1_delta_all[hidenode][hidenode], wh_2_delta_all[hidenode][hidenode], wh_3_delta_all[hidenode][hidenode], wh_4_delta_all[hidenode][hidenode];
		static double b_0_delta_all[hidenode], b_1_delta_all[hidenode], b_2_delta_all[hidenode], b_3_delta_all[hidenode], b_4_delta_all[hidenode];
		static double w_0_delta_all[col][hidenode], w_1_delta_all[hidenode][hidenode], w_2_delta_all[hidenode][hidenode], w_3_delta_all[hidenode][hidenode], w_4_delta_all[hidenode][hidenode];
		static double lamuda = atof(argv[2]);

		for (int epoch = 0; epoch < 5; epoch++)
		{
			srand(time(NULL));
			if (epoch != 0)
			{
				gaussrand();
			}
			for (int i = 0; i < hidenode; i++)
			{
				b_0[i] = 0;
				b_1[i] = 0;
				b_2[i] = 0;
				b_3[i] = 0;
				b_4[i] = 0;
				for (int j = 0; j < col; j++)
				{
					w_0[j][i] = gaussrand() * 0.2357;
				}
				for (int j = 0; j < hidenode; j++)
				{
					w_1[j][i] = gaussrand() * 0.254;
					w_2[j][i] = gaussrand() * 0.254;
					w_3[j][i] = gaussrand() * 0.254;
					w_4[j][i] = gaussrand() * 0.254;
				}
			}
			for (int i = 0; i < hidenode; i++)
			{
				for (int j = 0; j < hidenode; j++)
				{
					wh_0[i][j] = gaussrand() * 0.254;
					wh_1[i][j] = gaussrand() * 0.254;
					wh_2[i][j] = gaussrand() * 0.254;
					wh_3[i][j] = gaussrand() * 0.254;
					wh_4[i][j] = gaussrand() * 0.254;
				}
			}
			b1 = 0;
			for (int i = 0; i < hidenode; i++)
			{
				w1[i] = gaussrand() * 0.25;
			}
			//training start
			for (int i = 0; i < iter; i++)
			{
				alpha = 0.0005 * pow(0.5, (i / row) / 100.0);
				double y = outcome[i % row];
				//1,3,5,7,...,47
				for (int j = 0; j < timestep / 2; j++)
				{
					for (int k = 0; k < col; k++)
					{
						layer_0[k] = data[i % row][k][2 * j];
					}
					for (int k = 0; k < hidenode; k++)
					{
						double o1 = b_0[k];
						for (int m = 0; m < col; m++)
						{
							o1 += layer_0[m] * w_0[m][k];
						}
						o1_0[2 * j][k] = o1;
						layer_1_0[2 * j][k] = activation(o1);
					}
				}
				//2,4,6,8,...,48
				for (int j = 0; j < timestep / 2; j++)
				{
					for (int k = 0; k < col; k++)
					{
						layer_0[k] = data[i % row][k][2 * j + 1];
					}
					for (int k = 0; k < hidenode; k++)
					{
						double o1 = b_0[k];
						for (int m = 0; m < col; m++)
						{
							o1 += layer_0[m] * w_0[m][k];
						}
						for (int m = 0; m < hidenode; m++)
						{
							o1 += layer_1_0[2 * j][m] * wh_0[m][k];
						}
						o1_0[2 * j + 1][k] = o1;
						layer_1_0[2 * j + 1][k] = activation(o1);
					}
				}
				//2,6,10,14,...,46
				for (int j = 0; j < timestep / 4; j++)
				{
					for (int k = 0; k < hidenode; k++)
					{
						double o1 = b_1[k];
						for (int m = 0; m < hidenode; m++)
						{
							o1 += layer_1_0[4 * j + 1][m] * w_1[m][k];
						}
						o1_1[2 * j][k] = o1;
						layer_1_1[2 * j][k] = activation(o1);
					}
				}
				//4,8,12,16,...,48
				for (int j = 0; j < timestep / 4; j++)
				{
					for (int k = 0; k < hidenode; k++)
					{
						double o1 = b_1[k];
						for (int m = 0; m < hidenode; m++)
						{
							o1 += layer_1_0[4 * j + 3][m] * w_1[m][k];
						}
						for (int m = 0; m < hidenode; m++)
						{
							o1 += layer_1_1[2 * j][m] * wh_1[m][k];
						}
						o1_1[2 * j + 1][k] = o1;
						layer_1_1[2 * j + 1][k] = activation(o1);
					}
				}
				//4,12,20,28,36,44
				for (int j = 0; j < timestep / 8; j++)
				{
					for (int k = 0; k < hidenode; k++)
					{
						double o1 = b_2[k];
						for (int m = 0; m < hidenode; m++)
						{
							o1 += layer_1_1[4 * j + 1][m] * w_2[m][k];
						}
						o1_2[2 * j][k] = o1;
						layer_1_2[2 * j][k] = activation(o1);
					}
				}
				//8,16,24,32,40,48
				for (int j = 0; j < timestep / 8; j++)
				{
					for (int k = 0; k < hidenode; k++)
					{
						double o1 = b_2[k];
						for (int m = 0; m < hidenode; m++)
						{
							o1 += layer_1_1[4 * j + 3][m] * w_2[m][k];
						}
						for (int m = 0; m < hidenode; m++)
						{
							o1 += layer_1_2[2 * j][m] * wh_2[m][k];
						}
						o1_2[2 * j + 1][k] = o1;
						layer_1_2[2 * j + 1][k] = activation(o1);
					}
				}
				//8,24,40
				for (int j = 0; j < timestep / 16; j++)
				{
					for (int k = 0; k < hidenode; k++)
					{
						double o1 = b_3[k];
						for (int m = 0; m < hidenode; m++)
						{
							o1 += layer_1_2[4 * j + 1][m] * w_3[m][k];
						}
						o1_3[2 * j][k] = o1;
						layer_1_3[2 * j][k] = activation(o1);
					}
				}
				//16,32,48
				for (int j = 0; j < timestep / 16; j++)
				{
					for (int k = 0; k < hidenode; k++)
					{
						double o1 = b_3[k];
						for (int m = 0; m < hidenode; m++)
						{
							o1 += layer_1_2[4 * j + 3][m] * w_3[m][k];
						}
						for (int m = 0; m < hidenode; m++)
						{
							o1 += layer_1_3[2 * j][m] * wh_3[m][k];
						}
						o1_3[2 * j + 1][k] = o1;
						layer_1_3[2 * j + 1][k] = activation(o1);
					}
				}
				//16
				for (int j = 0; j < hidenode; j++)
				{
					double o1 = b_4[j];
					for (int k = 0; k < hidenode; k++)
					{
						o1 += layer_1_3[1][k] * w_4[k][j];
					}
					o1_4[0][j] = o1;
					layer_1_4[0][j] = activation(o1);
				}
				//32
				for (int j = 0; j < hidenode; j++)
				{
					double o1 = b_4[j];
					for (int k = 0; k < hidenode; k++)
					{
						o1 += layer_1_3[3][k] * w_4[k][j];
					}
					for (int k = 0; k < hidenode; k++)
					{
						o1 += layer_1_4[0][k] * wh_4[k][j];
					}
					o1_4[1][j] = o1;
					layer_1_4[1][j] = activation(o1);
				}
				//48
				for (int j = 0; j < hidenode; j++)
				{
					double o1 = b_4[j];
					for (int k = 0; k < hidenode; k++)
					{
						o1 += layer_1_3[5][k] * w_4[k][j];
					}
					for (int k = 0; k < hidenode; k++)
					{
						o1 += layer_1_4[1][k] * wh_4[k][j];
					}
					o1_4[2][j] = o1;
					layer_1_4[2][j] = activation(o1);
				}
				//Output layer
				double o2 = b1;
				for (int j = 0; j < hidenode; j++)
				{
					o2 += layer_1_4[2][j] * w1[j];
				}
				layer_2 = outactive(o2);
				layer_2_delta = (y / layer_2 - (1 - y) / (1 - layer_2)) * doutactive(layer_2);
				b1_delta = alpha * layer_2_delta;
				for (int j = 0; j < hidenode; j++)
				{
					w1_delta[j] = alpha * (layer_2_delta * layer_1_4[2][j] - lamuda * w1[j] / minibatchsize);
				}
				//48_4
				for (int j = 0; j < hidenode; j++)
				{
					layer_1_4_delta[2][j] = layer_2_delta * w1[j];
					layer_1_4_delta[2][j] = layer_1_4_delta[2][j] * dactivate(o1_4[2][j]);
					for (int k = 0; k < hidenode; k++)
					{
						wh_4_delta[k][j] = alpha * (layer_1_4_delta[2][j] * layer_1_4[1][k] - lamuda * wh_4[k][j] / minibatchsize);
					}
					b_4_delta[j] = alpha * layer_1_4_delta[2][j];
					for (int k = 0; k < hidenode; k++)
					{
						w_4_delta[k][j] = alpha * (layer_1_4_delta[2][j] * layer_1_3[5][k] - lamuda * w_4[k][j] / minibatchsize);
					}
				}
				//32_4
				for (int j = 0; j < hidenode; j++)
				{
					layer_1_4_delta[1][j] = 0.0;
					for (int k = 0; k < hidenode; k++)
					{
						layer_1_4_delta[1][j] += layer_1_4_delta[2][k] * wh_4[j][k];
					}
					layer_1_4_delta[1][j] = layer_1_4_delta[1][j] * dactivate(o1_4[1][j]);
					for (int k = 0; k < hidenode; k++)
					{
						wh_4_delta[k][j] += alpha * layer_1_4_delta[1][j] * layer_1_4[0][k];
					}
					b_4_delta[j] += alpha * layer_1_4_delta[1][j];
					for (int k = 0; k < hidenode; k++)
					{
						w_4_delta[k][j] += alpha * layer_1_4_delta[1][j] * layer_1_3[3][k];
					}
				}
				//16_4
				for (int j = 0; j < hidenode; j++)
				{
					layer_1_4_delta[0][j] = 0.0;
					for (int k = 0; k < hidenode; k++)
					{
						layer_1_4_delta[0][j] += layer_1_4_delta[1][k] * wh_4[j][k];
					}
					layer_1_4_delta[0][j] = layer_1_4_delta[0][j] * dactivate(o1_4[0][j]);
					b_4_delta[j] += alpha * layer_1_4_delta[0][j];
					for (int k = 0; k < hidenode; k++)
					{
						w_4_delta[k][j] += alpha * layer_1_4_delta[0][j] * layer_1_3[1][k];
					}
				}
				for (int j = 0; j < hidenode; j++)
				{
					b_3_delta[j] = 0;
					for (int k = 0; k < hidenode; k++)
					{
						wh_3_delta[j][k] = -alpha * lamuda * wh_3[j][k] / minibatchsize;
						w_3_delta[j][k] = -alpha * lamuda * w_3[j][k] / minibatchsize;
					}
				}
				//16,32,48
				for (int j = 0; j < timestep / 16; j++)
				{
					for (int k = 0; k < hidenode; k++)
					{
						layer_1_3_delta[2 * j + 1][k] = 0.0;
						for (int m = 0; m < hidenode; m++)
						{
							layer_1_3_delta[2 * j + 1][k] += layer_1_4_delta[j][m] * w_4[k][m];
						}
						layer_1_3_delta[2 * j + 1][k] = layer_1_3_delta[2 * j + 1][k] * dactivate(o1_3[2 * j + 1][k]);
						for (int m = 0; m < hidenode; m++)
						{
							wh_3_delta[m][k] += alpha * layer_1_3_delta[2 * j + 1][k] * layer_1_3[2 * j][m];
						}
						b_3_delta[k] += alpha * layer_1_3_delta[2 * j + 1][k];
						for (int m = 0; m < hidenode; m++)
						{
							w_3_delta[m][k] += alpha * layer_1_3_delta[2 * j + 1][k] * layer_1_2[4 * j + 3][m];
						}
					}
				}
				//8,24,40
				for (int j = 0; j < timestep / 16; j++)
				{
					for (int k = 0; k < hidenode; k++)
					{
						layer_1_3_delta[2 * j][k] = 0.0;
						for (int m = 0; m < hidenode; m++)
						{
							layer_1_3_delta[2 * j][k] += layer_1_3_delta[2 * j + 1][m] * wh_3[k][m];
						}
						layer_1_3_delta[2 * j][k] = layer_1_3_delta[2 * j][k] * dactivate(o1_3[2 * j][k]);
						b_3_delta[k] += alpha * layer_1_3_delta[2 * j][k];
						for (int m = 0; m < hidenode; m++)
						{
							w_3_delta[m][k] += alpha * layer_1_3_delta[2 * j][k] * layer_1_2[4 * j + 1][m];
						}
					}
				}
				for (int j = 0; j < hidenode; j++)
				{
					b_2_delta[j] = 0;
					for (int k = 0; k < hidenode; k++)
					{
						wh_2_delta[j][k] = -alpha * lamuda * wh_2[j][k] / minibatchsize;
						w_2_delta[j][k] = -alpha * lamuda * w_2[j][k] / minibatchsize;
					}
				}
				//8,16,24,32,40,48
				for (int j = 0; j < timestep / 8; j++)
				{
					for (int k = 0; k < hidenode; k++)
					{
						layer_1_2_delta[2 * j + 1][k] = 0.0;
						for (int m = 0; m < hidenode; m++)
						{
							layer_1_2_delta[2 * j + 1][k] += layer_1_3_delta[j][m] * w_3[k][m];
						}
						layer_1_2_delta[2 * j + 1][k] = layer_1_2_delta[2 * j + 1][k] * dactivate(o1_2[2 * j + 1][k]);
						for (int m = 0; m < hidenode; m++)
						{
							wh_2_delta[m][k] += alpha * layer_1_2_delta[2 * j + 1][k] * layer_1_2[2 * j][m];
						}
						b_2_delta[k] += alpha * layer_1_2_delta[2 * j + 1][k];
						for (int m = 0; m < hidenode; m++)
						{
							w_2_delta[m][k] += alpha * layer_1_2_delta[2 * j + 1][k] * layer_1_1[4 * j + 3][m];
						}
					}
				}
				//4,12,20,28,36,44
				for (int j = 0; j < timestep / 8; j++)
				{
					for (int k = 0; k < hidenode; k++)
					{
						layer_1_2_delta[2 * j][k] = 0.0;
						for (int m = 0; m < hidenode; m++)
						{
							layer_1_2_delta[2 * j][k] += layer_1_2_delta[2 * j + 1][m] * wh_2[k][m];
						}
						layer_1_2_delta[2 * j][k] = layer_1_2_delta[2 * j][k] * dactivate(o1_2[2 * j][k]);
						b_2_delta[k] += alpha * layer_1_2_delta[2 * j][k];
						for (int m = 0; m < hidenode; m++)
						{
							w_2_delta[m][k] += alpha * layer_1_2_delta[2 * j][k] * layer_1_1[4 * j + 1][m];
						}
					}
				}
				for (int j = 0; j < hidenode; j++)
				{
					b_1_delta[j] = 0;
					for (int k = 0; k < hidenode; k++)
					{
						wh_1_delta[j][k] = -alpha * lamuda * wh_1[j][k] / minibatchsize;
						w_1_delta[j][k] = -alpha * lamuda * w_1[j][k] / minibatchsize;
					}
				}
				//4,8,12,16,...,48
				for (int j = 0; j < timestep / 4; j++)
				{
					for (int k = 0; k < hidenode; k++)
					{
						layer_1_1_delta[2 * j + 1][k] = 0.0;
						for (int m = 0; m < hidenode; m++)
						{
							layer_1_1_delta[2 * j + 1][k] += layer_1_2_delta[j][m] * w_2[k][m];
						}
						layer_1_1_delta[2 * j + 1][k] = layer_1_1_delta[2 * j + 1][k] * dactivate(o1_1[2 * j + 1][k]);
						for (int m = 0; m < hidenode; m++)
						{
							wh_1_delta[m][k] += alpha * layer_1_1_delta[2 * j + 1][k] * layer_1_1[2 * j][m];
						}
						b_1_delta[k] += alpha * layer_1_1_delta[2 * j + 1][k];
						for (int m = 0; m < hidenode; m++)
						{
							w_1_delta[m][k] += alpha * layer_1_1_delta[2 * j + 1][k] * layer_1_0[4 * j + 3][m];
						}
					}
				}
				//2,6,10,14,...,46
				for (int j = 0; j < timestep / 4; j++)
				{
					for (int k = 0; k < hidenode; k++)
					{
						layer_1_1_delta[2 * j][k] = 0.0;
						for (int m = 0; m < hidenode; m++)
						{
							layer_1_1_delta[2 * j][k] += layer_1_1_delta[2 * j + 1][m] * wh_1[k][m];
						}
						layer_1_1_delta[2 * j][k] = layer_1_1_delta[2 * j][k] * dactivate(o1_1[2 * j][k]);
						b_1_delta[k] += alpha * layer_1_1_delta[2 * j][k];
						for (int m = 0; m < hidenode; m++)
						{
							w_1_delta[m][k] += alpha * layer_1_1_delta[2 * j][k] * layer_1_0[4 * j + 1][m];
						}
					}
				}
				for (int j = 0; j < hidenode; j++)
				{
					b_0_delta[j] = 0;
					for (int k = 0; k < hidenode; k++)
					{
						wh_0_delta[j][k] = -alpha * lamuda * wh_0[j][k] / minibatchsize;
					}
					for (int k = 0; k < col; k++)
					{
						w_0_delta[k][j] = -alpha * lamuda * w_0[k][j] / minibatchsize;
					}
				}
				//2,4,6,8,...,48
				for (int j = 0; j < timestep / 2; j++)
				{
					for (int k = 0; k < col; k++)
					{
						layer_0[k] = data[i % row][k][2 * j + 1];
					}
					for (int k = 0; k < hidenode; k++)
					{
						layer_1_0_delta[2 * j + 1][k] = 0.0;
						for (int m = 0; m < hidenode; m++)
						{
							layer_1_0_delta[2 * j + 1][k] += layer_1_1_delta[j][m] * w_1[k][m];
						}
						layer_1_0_delta[2 * j + 1][k] = layer_1_0_delta[2 * j + 1][k] * dactivate(o1_0[2 * j + 1][k]);
						for (int m = 0; m < hidenode; m++)
						{
							wh_0_delta[m][k] += alpha * layer_1_0_delta[2 * j + 1][k] * layer_1_0[2 * j][m];
						}
						b_0_delta[k] += alpha * layer_1_0_delta[2 * j + 1][k];
						for (int m = 0; m < col; m++)
						{
							w_0_delta[m][k] += alpha * layer_1_0_delta[2 * j + 1][k] * layer_0[m];
						}
					}
				}
				//1,3,5,7,...,47
				for (int j = 0; j < timestep / 2; j++)
				{
					for (int k = 0; k < col; k++)
					{
						layer_0[k] = data[i % row][k][2 * j];
					}
					for (int k = 0; k < hidenode; k++)
					{
						layer_1_0_delta[2 * j][k] = 0.0;
						for (int m = 0; m < hidenode; m++)
						{
							layer_1_0_delta[2 * j][k] += layer_1_0_delta[2 * j + 1][m] * wh_0[k][m];
						}
						layer_1_0_delta[2 * j][k] = layer_1_0_delta[2 * j][k] * dactivate(o1_0[2 * j][k]);
						b_0_delta[k] += alpha * layer_1_0_delta[2 * j][k];
						for (int m = 0; m < col; m++)
						{
							w_0_delta[m][k] += alpha * layer_1_0_delta[2 * j][k] * layer_0[m];
						}
					}
				}
				for (int j = 0; j < hidenode; j++)
				{
					b_0_delta_all[j] += b_0_delta[j] * (0.145475 + 0.854525 * y);
					b_1_delta_all[j] += b_1_delta[j] * (0.145475 + 0.854525 * y);
					b_2_delta_all[j] += b_2_delta[j] * (0.145475 + 0.854525 * y);
					b_3_delta_all[j] += b_3_delta[j] * (0.145475 + 0.854525 * y);
					b_4_delta_all[j] += b_4_delta[j] * (0.145475 + 0.854525 * y);
				}
				for (int j = 0; j < col; j++)
				{
					for (int k = 0; k < hidenode; k++)
					{
						w_0_delta_all[j][k] += w_0_delta[j][k] * (0.145475 + 0.854525 * y);
					}
				}
				for (int j = 0; j < hidenode; j++)
				{
					for (int k = 0; k < hidenode; k++)
					{
						w_1_delta_all[j][k] += w_1_delta[j][k] * (0.145475 + 0.854525 * y);
						w_2_delta_all[j][k] += w_2_delta[j][k] * (0.145475 + 0.854525 * y);
						w_3_delta_all[j][k] += w_3_delta[j][k] * (0.145475 + 0.854525 * y);
						w_4_delta_all[j][k] += w_4_delta[j][k] * (0.145475 + 0.854525 * y);
						wh_0_delta_all[j][k] += wh_0_delta[j][k] * (0.145475 + 0.854525 * y);
						wh_1_delta_all[j][k] += wh_1_delta[j][k] * (0.145475 + 0.854525 * y);
						wh_2_delta_all[j][k] += wh_2_delta[j][k] * (0.145475 + 0.854525 * y);
						wh_3_delta_all[j][k] += wh_3_delta[j][k] * (0.145475 + 0.854525 * y);
						wh_4_delta_all[j][k] += wh_4_delta[j][k] * (0.145475 + 0.854525 * y);
					}
				}
				b1_delta_all += b1_delta * (0.145475 + 0.854525 * y);
				for (int j = 0; j < hidenode; j++)
				{
					w1_delta_all[j] += w1_delta[j] * (0.145475 + 0.854525 * y);
				}
				if (i % minibatchsize == minibatchsize - 1)
				{
					for (int j = 0; j < hidenode; j++)
					{
						b_0[j] += b_0_delta_all[j];
						b_1[j] += b_1_delta_all[j];
						b_2[j] += b_2_delta_all[j];
						b_3[j] += b_3_delta_all[j];
						b_4[j] += b_4_delta_all[j];
						b_0_delta_all[j] = 0;
						b_1_delta_all[j] = 0;
						b_2_delta_all[j] = 0;
						b_3_delta_all[j] = 0;
						b_4_delta_all[j] = 0;
					}
					for (int j = 0; j < col; j++)
					{
						for (int k = 0; k < hidenode; k++)
						{
							w_0[j][k] += w_0_delta_all[j][k];
							w_0_delta_all[j][k] = 0;
						}
					}
					for (int j = 0; j < hidenode; j++)
					{
						for (int k = 0; k < hidenode; k++)
						{
							w_1[j][k] += w_1_delta_all[j][k];
							w_2[j][k] += w_2_delta_all[j][k];
							w_3[j][k] += w_3_delta_all[j][k];
							w_4[j][k] += w_4_delta_all[j][k];
							wh_0[j][k] += wh_0_delta_all[j][k];
							wh_1[j][k] += wh_1_delta_all[j][k];
							wh_2[j][k] += wh_2_delta_all[j][k];
							wh_3[j][k] += wh_3_delta_all[j][k];
							wh_4[j][k] += wh_4_delta_all[j][k];
							w_1_delta_all[j][k] = 0;
							w_2_delta_all[j][k] = 0;
							w_3_delta_all[j][k] = 0;
							w_4_delta_all[j][k] = 0;
							wh_0_delta_all[j][k] = 0;
							wh_1_delta_all[j][k] = 0;
							wh_2_delta_all[j][k] = 0;
							wh_3_delta_all[j][k] = 0;
							wh_4_delta_all[j][k] = 0;
						}
					}
					b1 += b1_delta_all;
					b1_delta_all = 0;
					for (int j = 0; j < hidenode; j++)
					{
						w1[j] += w1_delta_all[j];
						w1_delta_all[j] = 0;
					}
				}
				if (i % row == row - 1)
				{
					double e = 0;
					for (int index = 0; index < row; index++)
					{
						double y = outcome[index % row];
						//1,3,5,7,...,47
						for (int j = 0; j < timestep / 2; j++)
						{
							for (int k = 0; k < col; k++)
							{
								layer_0[k] = data[index % row][k][2 * j];
							}
							for (int k = 0; k < hidenode; k++)
							{
								double o1 = b_0[k];
								for (int m = 0; m < col; m++)
								{
									o1 += layer_0[m] * w_0[m][k];
								}
								layer_1_0[2 * j][k] = activation(o1);
							}
						}
						//2,4,6,8,...,48
						for (int j = 0; j < timestep / 2; j++)
						{
							for (int k = 0; k < col; k++)
							{
								layer_0[k] = data[index % row][k][2 * j + 1];
							}
							for (int k = 0; k < hidenode; k++)
							{
								double o1 = b_0[k];
								for (int m = 0; m < col; m++)
								{
									o1 += layer_0[m] * w_0[m][k];
								}
								for (int m = 0; m < hidenode; m++)
								{
									o1 += layer_1_0[2 * j][m] * wh_0[m][k];
								}
								layer_1_0[2 * j + 1][k] = activation(o1);
							}
						}
						//2,6,10,14,...,46
						for (int j = 0; j < timestep / 4; j++)
						{
							for (int k = 0; k < hidenode; k++)
							{
								double o1 = b_1[k];
								for (int m = 0; m < hidenode; m++)
								{
									o1 += layer_1_0[4 * j + 1][m] * w_1[m][k];
								}
								layer_1_1[2 * j][k] = activation(o1);
							}
						}
						//4,8,12,16,...,48
						for (int j = 0; j < timestep / 4; j++)
						{
							for (int k = 0; k < hidenode; k++)
							{
								double o1 = b_1[k];
								for (int m = 0; m < hidenode; m++)
								{
									o1 += layer_1_0[4 * j + 3][m] * w_1[m][k];
								}
								for (int m = 0; m < hidenode; m++)
								{
									o1 += layer_1_1[2 * j][m] * wh_1[m][k];
								}
								layer_1_1[2 * j + 1][k] = activation(o1);
							}
						}
						//4,12,20,28,36,44
						for (int j = 0; j < timestep / 8; j++)
						{
							for (int k = 0; k < hidenode; k++)
							{
								double o1 = b_2[k];
								for (int m = 0; m < hidenode; m++)
								{
									o1 += layer_1_1[4 * j + 1][m] * w_2[m][k];
								}
								layer_1_2[2 * j][k] = activation(o1);
							}
						}
						//8,16,24,32,40,48
						for (int j = 0; j < timestep / 8; j++)
						{
							for (int k = 0; k < hidenode; k++)
							{
								double o1 = b_2[k];
								for (int m = 0; m < hidenode; m++)
								{
									o1 += layer_1_1[4 * j + 3][m] * w_2[m][k];
								}
								for (int m = 0; m < hidenode; m++)
								{
									o1 += layer_1_2[2 * j][m] * wh_2[m][k];
								}
								layer_1_2[2 * j + 1][k] = activation(o1);
							}
						}
						//8,24,40
						for (int j = 0; j < timestep / 16; j++)
						{
							for (int k = 0; k < hidenode; k++)
							{
								double o1 = b_3[k];
								for (int m = 0; m < hidenode; m++)
								{
									o1 += layer_1_2[4 * j + 1][m] * w_3[m][k];
								}
								layer_1_3[2 * j][k] = activation(o1);
							}
						}
						//16,32,48
						for (int j = 0; j < timestep / 16; j++)
						{
							for (int k = 0; k < hidenode; k++)
							{
								double o1 = b_3[k];
								for (int m = 0; m < hidenode; m++)
								{
									o1 += layer_1_2[4 * j + 3][m] * w_3[m][k];
								}
								for (int m = 0; m < hidenode; m++)
								{
									o1 += layer_1_3[2 * j][m] * wh_3[m][k];
								}
								layer_1_3[2 * j + 1][k] = activation(o1);
							}
						}
						//16
						for (int j = 0; j < hidenode; j++)
						{
							double o1 = b_4[j];
							for (int k = 0; k < hidenode; k++)
							{
								o1 += layer_1_3[1][k] * w_4[k][j];
							}
							layer_1_4[0][j] = activation(o1);
						}
						//32
						for (int j = 0; j < hidenode; j++)
						{
							double o1 = b_4[j];
							for (int k = 0; k < hidenode; k++)
							{
								o1 += layer_1_3[3][k] * w_4[k][j];
							}
							for (int k = 0; k < hidenode; k++)
							{
								o1 += layer_1_4[0][k] * wh_4[k][j];
							}
							layer_1_4[1][j] = activation(o1);
						}
						//48
						for (int j = 0; j < hidenode; j++)
						{
							double o1 = b_4[j];
							for (int k = 0; k < hidenode; k++)
							{
								o1 += layer_1_3[5][k] * w_4[k][j];
							}
							for (int k = 0; k < hidenode; k++)
							{
								o1 += layer_1_4[1][k] * wh_4[k][j];
							}
							layer_1_4[2][j] = activation(o1);
						}
						double o2 = b1;
						for (int j = 0; j < hidenode; j++)
						{
							o2 += layer_1_4[2][j] * w1[j];
						}
						layer_2 = outactive(o2);
						e += -(y * log(layer_2) + (1 - y) * log(1 - layer_2)) * (0.145475 + 0.854525 * y);
					}
					outout << e;
					outout << ",";
					test(average, stdev, b_0, b_1, b_2, b_3, b_4, w_0, w_1, w_2, w_3, w_4, wh_0, wh_1, wh_2, wh_3, wh_4, b1, w1, result, 200 * epoch + i / (row * 1), loop, e);
					outout << e << endl;
					cout << e << endl;
				}
			}
		}
		ofstream out2("D:\\result_SRNN" + to_string(loop) + "_" + to_string(lamuda) + ".csv");
		for (int i = 0; i < 4000; i++)
		{
			for (int j = 0; j < 1000; j++)
			{
				out2 << result[i][j];
				out2 << ",";
			}
			out2 << endl;
		}
		out2.close();
		outout.close();
	}
	cout << "end" << endl;
}