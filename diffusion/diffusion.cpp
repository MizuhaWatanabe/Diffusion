// diffusion.cpp : �R���\�[�� �A�v���P�[�V�����̃G���g�� �|�C���g���`���܂��B
//

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N	  10000		// ���ԃX�e�b�v��
#define X	  100		// ���b�V����
#define time  1.0		// ���Ԃ̌v�Z�͈�
#define space 1.0		// ��Ԃ̌v�Z�͈�
#define D	  0.1		// �g�U�W��
#define write 100		// ���X�e�b�v�Ɉ�񏑂����ނ�

double u[N][X] = {};
double du[N][X] = {};
double x[X] = {};
double S[N][X] = {};
double mol[N] = {};
double t;
double dt, dx;
double a[4] = { 0, 0.5, 0.5 ,1 };
double b[4] = { (double)1 / 6, (double)1 / 3, (double)1 / 3, (double)1 / 6 };
double valid;
int    i, n, k, m;
FILE   *fp;
char   fname[32];

int _tmain(int argc, _TCHAR* argv[])
{
	// ���l���萫�̌v�Z�A�s����Ȃ�G���[
	dt = (double)time / N;
	dx = (double)space / X;
	valid = D * dt / (dx * dx);
	if (valid > 0.5){
		printf("ERROR : Calculation is unstable.(%e)\n", valid);
		getchar();
		return -1;
	}
	else{
		printf("Calculation is stable.(%e)\n", valid);
	}

	// �A�E�g�v�b�g�t�H���_�̏�����
	system("rd /s /q output");
	system("mkdir output");

	// ��������
	for (i = 0; i < X; i++){
		u[0][i] = 1;
	}

	// Runge-Kutta�@�̌v�Z
	for (n = 0; n < N; n++){

		//��Ԃ̏�����
		for (i = 0; i < X; i++){
			x[i] = i * dx;
		}

		// �A�E�g�v�b�g�t�@�C���̏����o��
		if (n % write == 0){
			sprintf(fname, "output/output_%d.csv", m);
			fp = fopen(fname, "w");
			if (fp == NULL){
				printf("cannot open output file.");
				return -1;
			}
			fprintf(fp, "x,u\n");
			for (i = 0; i < X; i++){
				fprintf(fp, "%e,%e\n", x[i], u[n][i]);
			}
			fclose(fp);
			m += 1;
		}
		
		// ���틫�E����
		u[n][X - 1] = exp(-5 * (double)n / N);

		// 4���ߎ��̌v�Z
		for (k = 0; k < 4; k++){
			for (i = 1; i < X - 1; i++){
				u[n][i] += a[k] * du[n][i];
				x[i] += a[k] * dx;
				du[n][i] = dt * D * ((u[n][i + 1] + u[n][i - 1] - 2 * u[n][i]) / (dx * dx) + (u[n][i + 1] - u[n][i - 1]) / (x[i] * dx));
				S[n][i] += b[k] * du[n][i];
			}
		}

		// ���틫�E����
		u[n][0] = 2 * u[n][1] - u[n][2];
		
		// ���Ԕ��W
		for (i = 0; i < X; i++){
			u[n + 1][i] = u[n][i] + S[n][i];
		}
	}

	//��Ԃ̏�����
	for (i = 0; i < X; i++){
		x[i] = i * dx;
	}

	// �S�����ʂ̎��ԕω��̌v�Z�B�Ō�̃Z���̓o���N��
	sprintf(fname, "output/mol.csv");
	fp = fopen(fname, "w");
	fprintf(fp, "t,totalMol\n");
	for (n = 0; n < N; n++){
		t = (double)n * dt;
		for (i = 0; i < X - 1; i++){
			mol[n] += 4 * M_PI * x[i] * x[i] *  u[n][i] * dx;
		}
		fprintf(fp, "%e,%e\n", t, mol[n]);
	}
	fclose(fp);

	// �v�Z�I���̊m�F
	printf("Calculation is complete. Please push enter key.");
	getchar();

	return 0;
}

