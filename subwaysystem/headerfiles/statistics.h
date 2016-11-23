#pragma once
//•W–{••Õ•ªU
double getmyu(double data[], int datasize) {
	double ans = 0;
	for (int i = 0; i < datasize;i++) {
		ans += data[i];
	}
	ans /= (datasize-1);
	return ans;
}

double getsoundmyu(double data[], int datasize) {
	double ans = 0;
	for (int i = 0; i < datasize; i++) {
		ans += fabs(data[i]);
	}
	ans /= (datasize-1);
	return ans;
}

double getVar(double data[], double myu, int datasize) {
	double ans = 0;
	for (int i = 0; i < datasize; i++) {
		ans += (data[i] - myu)* (data[i] - myu);
	}
	ans /= datasize;
	return ans;
}

double getSD(double data[], double myu, int datasize) {
	double ans = getVar(data, myu, datasize);
	ans = sqrt(ans);
	return ans;
}