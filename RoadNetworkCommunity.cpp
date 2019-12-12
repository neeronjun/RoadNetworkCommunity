// video_traj_analysis.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include<string>
#include<time.h>
#include <iostream>
#include<math.h>
#include<fstream>
#include<sstream>
#include <vector>
#include <stdlib.h>
#define pi 3.14159265358979323846

using namespace std;

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*::  This function converts decimal degrees to radians             :*/
/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
double deg2rad(double deg) {
	return (deg * pi / 180);
}

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*::  This function converts radians to decimal degrees             :*/
/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
double rad2deg(double rad) {
	return (rad * 180 / pi);
}

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*::  Function prototypes                                           :*/
/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
double deg2rad(double);
double rad2deg(double);

double distance(double lat1, double lon1, double lat2, double lon2, char unit) {
	double theta, dist;
	if ((lat1 == lat2) && (lon1 == lon2)) {
		return 0;
	}
	else {
		theta = lon1 - lon2;
		dist = sin(deg2rad(lat1)) * sin(deg2rad(lat2)) + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * cos(deg2rad(theta));
		dist = acos(dist);
		dist = rad2deg(dist);
		dist = dist * 60 * 1.1515;
		switch (unit) {
		case 'M':
			break;
		case 'K':
			dist = dist * 1.609344;
			break;
		case 'N':
			dist = dist * 0.8684;
			break;
		}
		return (dist);
	}
}

//min value cannot be zero.
void D2_calculation(int datasize, double space_length, int dimension, double min, double max, double interval) {

	clock_t time_req, time_receive;

	int iterations = static_cast<int>((max - min) / interval);
	int iteration_count = 0;
	//cout << "iterations: " << iterations << endl;
	double* fractal_dimensions = new double[2 * iterations];
	double sum_pi;

	for (int i = 0; i < 2 * iterations; i++)
		fractal_dimensions[i] = 0;

	double percentSum = 0;
	double D2 = 0;

	//char *inpname = "experimentationFile/synthetic_data/anti/R_10K.txt";
	//char *output = "experimentationFile/synthetic_data/anti/D2_value.txt";

	string inpname = "C:/Users/nrai/Desktop/samples30.txt";
	string output = "C:/Users/nrai/Desktop/D2_value.txt";

	double cell_side_r;
	double interval_num_f;
	int cell_count;
	double* cells = NULL;

	string line;
	string subline;

	int count;
	double a, b, c, d, e, f;

	double temp;
	//int indexTemp;
	int index_position;

	time_req = clock();

	if (dimension == 2) {
		for (double r = min; r <= max; r += interval) {
			sum_pi = 0;
			cell_side_r = r / space_length;
			interval_num_f = ceil(space_length / r);
			cell_count = (int)pow(interval_num_f, dimension);
			cells = new double[cell_count];

			std::ifstream in(inpname);

			for (int index = 0; index < cell_count; index++)
				cells[index] = 0;

			if (in)
			{
				while (getline(in, line))
				{
					count = 0;
					a = 0;
					b = 0;
					stringstream ss(line);
					while (ss >> subline)
					{
						if (count == 1)
						{
							temp = atof(subline.c_str());
							a = floor(temp / r);
							if (a >= interval_num_f) a = interval_num_f - 1;
						}
						if (count == 2)
						{
							temp = atof(subline.c_str());
							b = floor(temp / r);
							if (b >= interval_num_f) b = interval_num_f - 1;
						}

						count++;
					}
					index_position = (int)(a * pow(interval_num_f, 1) + b);
					cells[index_position]++;
				}
			}
			else
			{
				cout << "no such file" << endl;
			}

			in.close();

			for (int index = 0; index < cell_count; index++)
				sum_pi += pow(cells[index] / datasize, 2);

			fractal_dimensions[2 * iteration_count] = log(sum_pi);
			cout << log(sum_pi) << endl;
			fractal_dimensions[2 * iteration_count + 1] = log(cell_side_r);
			cout << log(cell_side_r) << endl;
			iteration_count++;

			//cout << "r = " << r << endl;
			cout << "count = " << iteration_count << endl;

			delete[] cells;
			cells = NULL;
		}

		/*
		int comparison_count = 0;
		for (int i = 0; i < iterations-1; i++)
		for (int j = i + 1; j < iterations;j++)

		{
		comparison_count++;
		D2 += (fractal_dimensions[2  j] - fractal_dimensions[2  i]) / ((fractal_dimensions[2  j + 1] - fractal_dimensions[2  i + 1]));
		}
		*/

		/*
		int comparison_count = 0;
		for (int i = 0; i < 25; i++){
		D2 += (fractal_dimensions[2  i] - fractal_dimensions[0]) / ((fractal_dimensions[2  i + 1] - fractal_dimensions[1]));
		comparison_count;
		}
		*/


		//D2 /= comparison_count;



		//D2 = abs((fractal_dimensions[2 * 6] - fractal_dimensions[2 * 7]) / ((fractal_dimensions[2 * 6 + 1] - fractal_dimensions[2 * 7 + 1])));
		D2 = abs((fractal_dimensions[2 * 0] - fractal_dimensions[2 * 1]) / ((fractal_dimensions[2 * 0 + 1] - fractal_dimensions[2 * 1 + 1])));
		std::cout << "D2 = " << D2 << endl;

		ofstream fout(output);

		//fout << "comparison_count= " << comparison_count << "D2= " << D2 << endl;

		for (int i = 0; i < iterations; i++)
			fout << fractal_dimensions[2 * i] << " : " << fractal_dimensions[2 * i + 1] << endl;

		fout << "D2= " << D2 << endl;

		fout.close();

		time_receive = clock();

		std::cout << "it took " << (time_receive - time_req) << "clock to run" << endl;

	}

	//delete [] fractal_dimensions;
	//fractal_dimensions = NULL;

}

void estimate_neighbor() {
	float nb, D2, r, shape[2], square[2], vol_shape, vol_sq;
	int N = 30000;

	vol_shape = 1;
	vol_sq = 2;
	r = 0.1;
	D2 = 1.99765;

	nb = powf((vol_shape / vol_sq), D2 / 2) * (N - 1) * powf(2, D2) * powf(r, D2);
	cout << nb;
}

void read_csv(double x1[], double y1[], double x2[], double y2[]) {
	fstream fin;
	fin.open("C:/Users/nrai/Downloads/trips1.csv", ios::in);
	double x = 0;
	int count = 0, count1 = 0;
	string trip_id;
	float segment[30][2];
	//segment[5][0] = { 0};

	//vector<string> row;
	string line, word, temp;

	while (fin >> temp) {
		//row.clear();
		getline(fin, line);

		stringstream s(line);
		count1 = 0;
		while (getline(s, word, ',')) {
			if (count1 == 0) {
				trip_id = word;
			}
			else if (count1 > 0 && count1 < 10) {
				x = atof(word.c_str());
				//cout << word << "," << x << endl;
			}
			else {
				break;
			}
			//cout << word << endl;
			count1++;
		}

		//trip_id = stoi(row[0]);
		//cout << trip_id << endl;
		count++;
		if (count == 3) {
			break;
		}
	}
	fin.close();

}

double distance(double x1, double y1, double x2, double y2) {
	double d = 0;
	d = sqrt((pow((x2 - x1), 2) + pow((y2 - y1), 2)));
	return d * 1000;
}

void generate_samples(double start_x, double start_y, double end_x, double end_y) {

	double temp_x, temp_y, inc, t, dt, dist = 0, count = 0, a_x, a_y;
	ofstream outfile;
	outfile.open("C:/Users/nrai/Desktop/samples30.txt", ios::in | ios::app); // append instead of overwrite

	dist = distance(start_x, start_y, end_x, end_y);
	dt = dist / 10;
	inc = dt;
	if (dist < dt) {
		cout << "small distance: " << dist << endl;
	}
	else {
		while (dt < dist) {
			t = dt / dist;
			//cout << "dt: " << dt << endl;
			//t = dt / d;
			temp_x = ((1 - t) * start_x + t * end_x);
			temp_y = ((1 - t) * start_y + t * end_y);
			dt += inc;
			outfile.precision(18);
			//if (temp_x > 40.9 && temp_x < 40.9 && temp_y > -74 && temp_y < -72.9)
			//	outfile << count++ <<", " << temp_x << ", " << temp_y << endl;
			a_x = 10 * (temp_x - 29.21) / (43.75 - 29.21) + 0;
			a_y = 10 * (temp_y + 122.73) / (-71.9 + 122.73) + 0;
			//if (a_x > 0 && a_x < 10 && a_y > 0 && a_y < 10)
			outfile << count++ << ", " << a_x << ", " << a_y << endl;
			//cout << count << ", " << a_x << ", " << a_y << endl;
			//cout << "SDFSD:" << distance(start_x, start_y, temp_x, temp_y) << endl;

		}
		//outfile << "asdfasd" << endl;
	}
	outfile.close();
}

int main()
{
	int datasize = 10000, dimension = 2, count1 = 0, count = 0, j = 0, n = 0;
	double min = 0.1, max = 10, interval = 1, space_length = 10, length = 3;
	double dist, x1, x2, y1, y2;
	double* traj = new double[500];
	double start_x, start_y, end_x, end_y;
	double center_x, center_y, sum_dist = 0, dt, t, a_x, a_y, len_max = 1;
	fstream fin, fcen;
	int i = 0, x = 0;
	fin.open("C:/Users/nrai/Desktop/trips-new.txt", ios::in);
	fcen.precision(18);
	fcen.open("C:/Users/nrai/Desktop/samples30.txt", ios::in | ios::app);
	string line, word, temp;

	cout.precision(18); 
	/*while (fin >> temp) {
		j = 0;
		count = 0;
		getline(fin, line);

		stringstream s(line);
		//count1 = 0;
		while (getline(s, word, '\t')) {
			if (word == "")
				continue;
			traj[count] = atof(word.c_str());
			//cout << word << endl;
			count++;
		}
		n = (int) traj[j++];
		count1++;
		//j = j + 2;
		//cout << "n: " << n << endl;
		sum_dist = 0;
		//m = m + n;
		//cout << "m: " << m << endl;
		x1 = traj[j++], y1 = traj[j++];
		//cout << "x1:" << x1 << "y1:" << y1 << "x2:" << traj[n * 2 - 1] << "y2:" << traj[2 * n] << endl;
		//cout << distance(x1, y1, traj[2 * n - 1], traj[n * 2]) << endl;
		//generate_samples(x1, y1, traj[2 * n - 1], traj[2 * n]);

		start_x = x1, start_y = y1;
		//calculate centers
		for (int k = 0; k < n - 1; k++) {
			x2 = traj[j++];
			y2 = traj[j++];

			dist = distance(x1, y1, x2, y2);
			//cout << "distnace: " << dist << endl;
			sum_dist += dist;
			//cout << "sum_dist: " << sum_dist << endl;
			if (sum_dist > 0.5) {

				sum_dist = sum_dist - dist;
				dt = length - sum_dist;
				t = dt / dist;

				//cout << "dt: " << dt << endl;
				//t = dt / d;
				end_x = ((1 - t) * x1 + t * x2);
				end_y = ((1 - t) * y1 + t * y2);


				center_x = (start_x + end_x) / 2;
				center_y = (start_y + end_y) / 2;

				a_x = 10 * (center_x - 29.21) / (43.75 - 29.21) + 0;
				a_y = 10 * (center_y + 122.73) / (-71.9 + 122.73) + 0;
				//if (a_x > 0 && a_x < 10 && a_y > 0 && a_y < 10)
				fcen << count << ", " << a_x << ", " << a_y << endl;

				generate_samples(start_x, start_y, end_x, end_y);
				start_x = end_x, start_y = end_y;
				sum_dist = 0;
				x2 = end_x, y2 = end_y;
			}

			if (k == n - 2) {
				//cout << "n:::" << distance(start_x, start_y, x2, y2) << endl;
				center_x = (start_x + x2) / 2;
				center_y = (start_y + y2) / 2;
				a_x = 10 * (center_x - 29.21) / (43.75 - 29.21) + 0;
				a_y = 10 * (center_y + 122.73) / (-71.9 + 122.73) + 0;
				//if (a_x > 0 && a_x < 10 && a_y > 0 && a_y < 10)
				fcen << count << ", " << a_x << ", " << a_y << endl;
				//cout << "n::kk:" << distance(centers[0], centers[1], start_x, start_y) << endl;
			}
			x1 = x2, y1 = y2;

		}
		//cout << count << endl;
	}
	*/
	std::cout << "Hello World!\n";

	//estimate_neighbor();	

	D2_calculation(datasize, space_length, dimension, min, max, interval);
	//fcen.close();
	fin.close();
	delete[] traj;
	return 0;
}
