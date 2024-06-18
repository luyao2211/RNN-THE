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
	//return 0.0042 * x * x + 0.5 * x;
	return 0.5 * (x + abs(x));
	//return tanh(x);
}

double dactivate(double x)
{
	//return 0.0042 * 2 * x + 0.5;
	if (x < 0)
		return 0;
	else
		return 1;
	//return 1 - tanh(x) * tanh(x);
}

double outactive(double x)
{
	int e = rand() * 11.0 / RAND_MAX - 5.5;
	if ((x + e - round(e)) < -10)
		return 4.5397868702434395e-05;
	else if ((x + e - round(e)) < -9)
		return 7.7996707283797329e-05 * x + 0.00082536494154040769;
	else if ((x + e - round(e)) < -8)
		return 0.00021195555448024638 * x + 0.0020309945663084493;
	else if ((x + e - round(e)) < -7)
		return 0.00057570106393416728 * x + 0.0049409586419398160;
	else if ((x + e - round(e)) < -6)
		return 0.0015615719622341291 * x + 0.011842054930039550;
	else if ((x + e - round(e)) < -5)
		return 0.0042202277676500807 * x + 0.027793989762535259;
	else if ((x + e - round(e)) < -4)
		return 0.011293359037806703 * x + 0.063159646113318368;
	else if ((x + e - round(e)) < -3)
		return 0.029439663215475222 * x + 0.13574486282399245;
	else if ((x + e - round(e)) < -2)
		return 0.071777048844550773 * x + 0.26275701971121912;
	else if ((x + e - round(e)) < -1)
		return 0.14973849934787756 * x + 0.41867992071787263;
	else if ((x + e - round(e)) < 0)
		return 0.23105857863000490 * x + 0.5;
	else if ((x + e - round(e)) < 1)
		return 0.23105857863000490 * x + 0.5;
	else if ((x + e - round(e)) < 2)
		return 0.14973849934787742 * x + 0.58132007928212748;
	else if ((x + e - round(e)) < 3)
		return 0.071777048844551050 * x + 0.73724298028878021;
	else if ((x + e - round(e)) < 4)
		return 0.029439663215475087 * x + 0.86425513717600810;
	else if ((x + e - round(e)) < 5)
		return 0.011293359037806816 * x + 0.93684035388668119;
	else if ((x + e - round(e)) < 6)
		return 0.0042202277676500755 * x + 0.97220601023746489;
	else if ((x + e - round(e)) < 7)
		return 0.0015615719622340540 * x + 0.98815794506996102;
	else if ((x + e - round(e)) < 8)
		return 0.00057570106393423082 * x + 0.99505904135805978;
	else if ((x + e - round(e)) < 9)
		return 0.00021195555448005887 * x + 0.99796900543369316;
	else if ((x + e - round(e)) < 10)
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

void test(double average[col], double stdev[col], double b[hidenode], double w[col][hidenode], double wh[hidenode][hidenode], double b1, double w1[hidenode], double(&result)[4000][1000], int index, int loop, double& e)
{
	ifstream in3("D:\\TestTestData" + to_string(loop) + ".csv");
	static double testData[4000][col][timestep];
	int d1 = 0, a1 = 0;
	string line;
	std::getline(in3, line);
	while (std::getline(in3, line))
	{
		for (int i = 0; i < timestep; i++)
		{
			for (int j = 0; j < col; j++)
			{
				a1 = (int)line.find(',');
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
	static double layer_1_pre[hidenode];
	static double layer_1_vector[timestep + 1][hidenode];
	double layer_2;
	//loss
	e = 0;
	for (int i = 0; i < 4000; i++)
	{
		for (int j = 0; j < timestep; j++)
		{
			for (int k = 0; k < col; k++)
			{
				layer_0[k] = testData[i][k][j];
			}
			for (int k = 0; k < hidenode; k++)
			{
				double o1 = b[k];
				for (int m = 0; m < col; m++)
				{
					o1 += layer_0[m] * w[m][k];
				}
				for (int m = 0; m < hidenode; m++)
				{
					layer_1_pre[m] = layer_1_vector[j][m];
				}
				for (int m = 0; m < hidenode; m++)
				{
					o1 += layer_1_pre[m] * wh[m][k];
				}
				layer_1_vector[j + 1][k] = activation(o1);
			}
		}
		//Output layer
		double o2 = b1;
		for (int j = 0; j < hidenode; j++)
		{
			o2 += layer_1_vector[timestep][j] * w1[j];
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
		ofstream outout("D:\\Convergence_RNN" + to_string(loop) + "_" + to_string(atof(argv[2])) + ".csv");
		//data preparation
		ifstream in1("D:\\TestData" + to_string(loop) + ".csv");
		static double data[row][col][timestep];
		string line;
		size_t a1, d1 = 0;
		std::getline(in1, line);
		while (std::getline(in1, line))
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
		std::getline(in2, line);
		for (int i = 0; i < row; i++)
		{
			std::getline(in2, line);
			outcome[i] = atof(line.c_str());
		}
		in2.close();
		static double result[4000][1000];
		//training preparation
		int iter = row * 200;
		int minibatchsize = 250;
		double alpha;
		static double b[hidenode];
		static double w[col][hidenode];
		static double wh[hidenode][hidenode];
		static double b1;
		static double w1[hidenode];
		static double layer_0[col];
		static double o1[timestep][hidenode];
		static double layer_1_vector[timestep + 1][hidenode];
		static double layer_1[hidenode];
		static double layer_1_pre[hidenode];
		double layer_2;
		double layer_2_delta;
		static double layer_1_delta[hidenode];
		static double layer_1_future_delta[hidenode];
		double b1_delta;
		static double w1_delta[hidenode];
		static double wh_delta[hidenode][hidenode];
		static double b_delta[hidenode];
		static double w_delta[col][hidenode];
		double b1_delta_all = 0;
		static double w1_delta_all[hidenode];
		static double wh_delta_all[hidenode][hidenode];
		static double b_delta_all[hidenode];
		static double w_delta_all[col][hidenode];
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
				b[i] = 0;
				for (int j = 0; j < col; j++)
				{
					w[j][i] = gaussrand() * 0.17277;
				}
			}
			for (int i = 0; i < hidenode; i++)
			{
				for (int j = 0; j < hidenode; j++)
				{
					wh[i][j] = gaussrand() * 0.1796;
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
				for (int j = 0; j < timestep; j++)
				{
					for (int k = 0; k < col; k++)
					{
						layer_0[k] = data[i % row][k][j];
					}
					for (int k = 0; k < hidenode; k++)
					{
						o1[j][k] = b[k];
						for (int m = 0; m < col; m++)
						{
							o1[j][k] += layer_0[m] * w[m][k];
						}
						for (int m = 0; m < hidenode; m++)
						{
							layer_1_pre[m] = layer_1_vector[j][m];
						}
						for (int m = 0; m < hidenode; m++)
						{
							o1[j][k] += layer_1_pre[m] * wh[m][k];
						}
						layer_1_vector[j + 1][k] = activation(o1[j][k]);
					}
				}
				//Output layer
				double o2 = b1;
				for (int j = 0; j < hidenode; j++)
				{
					o2 += layer_1_vector[timestep][j] * w1[j];
				}
				layer_2 = outactive(o2);
				layer_2_delta = (y / layer_2 - (1 - y) / (1 - layer_2)) * doutactive(layer_2);//(y - layer_2) * doutactive(layer_2);
				b1_delta = alpha * layer_2_delta;
				for (int j = 0; j < hidenode; j++)
				{
					w1_delta[j] = alpha * (layer_2_delta * layer_1_vector[timestep][j] - lamuda * w1[j] / minibatchsize);
				}
				for (int j = 0; j < hidenode; j++)
				{
					layer_1_delta[j] = layer_2_delta * w1[j];
					layer_1_delta[j] = layer_1_delta[j] * dactivate(o1[timestep - 1][j]);
					for (int k = 0; k < hidenode; k++)
					{
						wh_delta[k][j] = alpha * (layer_1_delta[j] * layer_1_vector[timestep - 1][k] - lamuda * wh[k][j] / minibatchsize);
					}
					b_delta[j] = alpha * layer_1_delta[j];
					for (int k = 0; k < col; k++)
					{
						w_delta[k][j] = alpha * (layer_1_delta[j] * layer_0[k] - lamuda * w[k][j] / minibatchsize);
					}
				}
				for (int j = 0; j < hidenode; j++)
				{
					layer_1_future_delta[j] = layer_1_delta[j];
				}
				for (int j = timestep - 2; j > -1; j--)
				{
					for (int k = 0; k < col; k++)
					{
						layer_0[k] = data[i % row][k][j];
					}
					for (int k = 0; k < hidenode; k++)
					{
						layer_1_pre[k] = layer_1_vector[j][k];
					}
					for (int k = 0; k < hidenode; k++)
					{
						layer_1_delta[k] = 0.0;
						for (int m = 0; m < hidenode; m++)
						{
							layer_1_delta[k] += layer_1_future_delta[m] * wh[k][m];
						}
						layer_1_delta[k] = layer_1_delta[k] * dactivate(o1[j][k]);
						for (int m = 0; m < hidenode; m++)
						{
							wh_delta[m][k] += alpha * layer_1_delta[k] * layer_1_pre[m];
						}
						b_delta[k] += alpha * layer_1_delta[k];
						for (int m = 0; m < col; m++)
						{
							w_delta[m][k] += alpha * layer_1_delta[k] * layer_0[m];
						}
					}
					for (int k = 0; k < hidenode; k++)
					{
						layer_1_future_delta[k] = layer_1_delta[k];
					}
				}
				for (int j = 0; j < hidenode; j++)
				{
					b_delta_all[j] += b_delta[j] * (0.145475 + 0.854525 * y);
				}
				for (int j = 0; j < col; j++)
				{
					for (int k = 0; k < hidenode; k++)
					{
						w_delta_all[j][k] += w_delta[j][k] * (0.145475 + 0.854525 * y);
					}
				}
				for (int j = 0; j < hidenode; j++)
				{
					for (int k = 0; k < hidenode; k++)
					{
						wh_delta_all[j][k] += wh_delta[j][k] * (0.145475 + 0.854525 * y);
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
						b[j] += b_delta_all[j];
						b_delta_all[j] = 0;
					}
					for (int j = 0; j < col; j++)
					{
						for (int k = 0; k < hidenode; k++)
						{
							w[j][k] += w_delta_all[j][k];
							w_delta_all[j][k] = 0;
						}
					}
					for (int j = 0; j < hidenode; j++)
					{
						for (int k = 0; k < hidenode; k++)
						{
							wh[j][k] += wh_delta_all[j][k];
							wh_delta_all[j][k] = 0;
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
					//loss
					double e = 0;
					for (int index = 0; index < row; index++)
					{
						double y = outcome[index];
						for (int j = 0; j < timestep; j++)
						{
							for (int k = 0; k < col; k++)
							{
								layer_0[k] = data[index][k][j];
							}
							for (int k = 0; k < hidenode; k++)
							{
								double o1 = b[k];
								for (int m = 0; m < col; m++)
								{
									o1 += layer_0[m] * w[m][k];
								}
								for (int m = 0; m < hidenode; m++)
								{
									layer_1_pre[m] = layer_1_vector[j][m];
								}
								for (int m = 0; m < hidenode; m++)
								{
									o1 += layer_1_pre[m] * wh[m][k];
								}
								layer_1_vector[j + 1][k] = activation(o1);
							}
						}
						double o2 = b1;
						for (int j = 0; j < hidenode; j++)
						{
							o2 += layer_1_vector[timestep][j] * w1[j];
						}
						layer_2 = outactive(o2);
						e += -(y * log(layer_2) + (1 - y) * log(1 - layer_2)) * (0.145475 + 0.854525 * y);
					}
					outout << e;
					outout << ",";
					test(average, stdev, b, w, wh, b1, w1, result, 200 * epoch + i / (row * 1), loop, e);
					outout << e << endl;
					cout << e << endl;
				}
			}
		}
		ofstream out2("D:\\result_RNN" + to_string(loop) + "_" + to_string(lamuda) + ".csv");
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