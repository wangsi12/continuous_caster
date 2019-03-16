#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define Section 12
#define NNN 4001

void Physicial_Parameters(float);
float Boundary_Condition(int j,int ny,float Ly,float *,float *);
void Monitor_MeanTemperature(int , int , int , int *);
void OutputTemperature(int nx,int ny,int tnpts,int out_num);

float T_Last[11][NNN];
float T_New[11][NNN];
float T[11][NNN][21] = { 0 };
float Vcast=-0.02,h=1000.0,lamda=50.0,Ce=540.0,pho=7000.0,a=1.0;

int main()
{
	int nx=11,ny=NNN,tnpts=10001,out_num=500,count=0;;
	float Lx=0.125,Ly=28.599,t_final=2000.0,T_Cast=1558,Tw=30, dx,dy,tao;
	float ccml[Section+1]={0.0,0.2,0.4,0.6,0.8,1.0925,2.27,4.29,5.831,9.6065,13.6090,19.87014,28.599};
	float H_Init[Section] = { 1380,1170,980,800,1223.16,735.05,424.32,392.83,328.94,281.64,246.16,160.96 };
	float y[30001];
	int Y_Label[Section+1]={0};
	float T_Up,T_Down,T_Left,T_Right;
	bool flag = true;
	
	clock_t start,end,AssumeTime;

	dx=Lx/(nx-1);
	dy=Ly/(ny-1);
	tao=t_final/(tnpts-1);

	for(int j=0;j<ny;j++)
	{
		y[j]=j*dy;
	}

	Y_Label[Section]=ny-1;
	for(int j=0;j<ny-1;j++)
	{
		for(int i=1;i<Section;i++)
		{
			if(y[j]<=ccml[i]&&y[j+1]>=ccml[i])
				Y_Label[i]=j+1;
		}
	}

	printf("Casting Temperature = %f ", T_Cast);
	printf("\n");
	printf("The thick of steel billets(m) = %f ", Lx);
	printf("\n");
	printf("The length of steel billets(m) = %f ", Ly);
	printf("\n");
	printf("dx(m) = %f ", dx);
	printf("dy(m) = %f ", dy);
	printf("tao(s) = %f ", tao);
	printf("\n");
	printf("simulation time(s) = %f ", t_final);

	start=clock();//开始计时
	for(int i=0;i<nx;i++)
	{
		for(int j=0;j<ny;j++)
		{
			T_Last[i][j]=T_Cast;//初始化温度
		}
	}
	
	for(int k=0;k<tnpts-1;k++)
	{
		for(int j=0;j<ny;j++)
		{
			//h=Boundary_Condition(j,ny,Ly,ccml,H_Init);			
			for(int i=0;i<nx;i++)
			{
				if (flag)
				{
					Physicial_Parameters(T_Last[i][j]);
					a = lamda / (pho*Ce);
					if (i == 0 && j != 0 && j != ny - 1)//散热边，下边
					{
						h = Boundary_Condition(j, ny, Ly, ccml, H_Init);
						T_Up = T_Last[i + 1][j];
						T_Down = T_Last[i + 1][j] - 2 * dx*h*(T_Last[i][j] - Tw) / lamda;
						T_Right = T_Last[i][j + 1];
						T_Left = T_Last[i][j - 1];
						T_New[i][j] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + (1 - 2 * a*(tao / (dx*dx)) - 2 * a*(tao / (dy*dy)) + tao*Vcast / dy)*T_Last[i][j] +
							a*(tao / (dy*dy))*T_Right + (a*(tao / (dy*dy)) - tao*Vcast / dy)*T_Left;
					}
					else if (i == nx - 1 && j != 0 && j != ny - 1)//上边
					{
						T_Up = T_Last[i - 1][j];
						T_Down = T_Last[i - 1][j];
						T_Right = T_Last[i][j + 1];
						T_Left = T_Last[i][j - 1];
						T_New[i][j] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + (1 - 2 * a*(tao / (dx*dx)) - 2 * a*(tao / (dy*dy)) + tao*Vcast / dy)*T_Last[i][j] +
							a*(tao / (dy*dy))*T_Right + (a*(tao / (dy*dy)) - tao*Vcast / dy)*T_Left;
					}
					else if (j == 0 && i != 0 && i != nx - 1)//左边
					{
						T_New[i][j] = T_Cast;
					}
					else if (j == ny - 1 && i != 0 && i != nx - 1)//右边
					{
						T_Up = T_Last[i + 1][j];
						T_Down = T_Last[i - 1][j];
						T_Right = T_Last[i][j - 1];
						T_Left = T_Last[i][j - 1];
						T_New[i][j] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + (1 - 2 * a*(tao / (dx*dx)) - 2 * a*(tao / (dy*dy)) + tao*Vcast / dy)*T_Last[i][j] +
							a*(tao / (dy*dy))*T_Right + (a*(tao / (dy*dy)) - tao*Vcast / dy)*T_Left;
					}
					else if (i == 0 && j == 0)//左下角
					{
						T_New[i][j] = T_Cast;
					}
					else if (i == 0 && j == ny - 1)//右下角
					{
						h = Boundary_Condition(j, ny, Ly, ccml, H_Init);
						T_Up = T_Last[i + 1][j];
						T_Down = T_Last[i + 1][j] - 2 * dx*h*(T_Last[i][j] - Tw) / lamda;
						T_Right = T_Last[i + 1][j];
						T_Left = T_Last[i][j - 1];
						T_New[i][j] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + (1 - 2 * a*(tao / (dx*dx)) - 2 * a*(tao / (dy*dy)) + tao*Vcast / dy)*T_Last[i][j] +
							a*(tao / (dy*dy))*T_Right + (a*(tao / (dy*dy)) - tao*Vcast / dy)*T_Left;
					}
					else if (i == nx - 1 && j == 0)//左上角
					{
						T_New[i][j] = T_Cast;
					}
					else if (i == nx - 1 && j == ny - 1)//右上角
					{
						T_Up = T_Last[i - 1][j];
						T_Down = T_Last[i - 1][j];
						T_Right = T_Last[i - 1][j];
						T_Left = T_Last[i][j - 1];
						T_New[i][j] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + (1 - 2 * a*(tao / (dx*dx)) - 2 * a*(tao / (dy*dy)) + tao*Vcast / dy)*T_Last[i][j] +
							a*(tao / (dy*dy))*T_Right + (a*(tao / (dy*dy)) - tao*Vcast / dy)*T_Left;
					}
					else//内部
					{
						T_Up = T_Last[i + 1][j];
						T_Down = T_Last[i - 1][j];
						T_Right = T_Last[i][j + 1];
						T_Left = T_Last[i][j - 1];
						T_New[i][j] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + (1 - 2 * a*(tao / (dx*dx)) - 2 * a*(tao / (dy*dy)) + tao*Vcast / dy)*T_Last[i][j] +
							a*(tao / (dy*dy))*T_Right + (a*(tao / (dy*dy)) - tao*Vcast / dy)*T_Left;
					}
				}
				else
				{
					Physicial_Parameters(T_New[i][j]);
					a = lamda / (pho*Ce);
					if (i == 0 && j != 0 && j != ny - 1)//散热边，下边
					{
						h = Boundary_Condition(j, ny, Ly, ccml, H_Init);
						T_Up = T_New[i + 1][j];
						T_Down = T_New[i + 1][j] - 2 * dx*h*(T_New[i][j] - Tw) / lamda;
						T_Right = T_New[i][j + 1];
						T_Left = T_New[i][j - 1];
						T_Last[i][j] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + (1 - 2 * a*(tao / (dx*dx)) - 2 * a*(tao / (dy*dy)) + tao*Vcast / dy)*T_Last[i][j] +
							a*(tao / (dy*dy))*T_Right + (a*(tao / (dy*dy)) - tao*Vcast / dy)*T_Left;
					}
					else if (i == nx - 1 && j != 0 && j != ny - 1)//上边
					{
						T_Up = T_New[i - 1][j];
						T_Down = T_New[i - 1][j];
						T_Right = T_New[i][j + 1];
						T_Left = T_New[i][j - 1];
						T_Last[i][j] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + (1 - 2 * a*(tao / (dx*dx)) - 2 * a*(tao / (dy*dy)) + tao*Vcast / dy)*T_Last[i][j] +
							a*(tao / (dy*dy))*T_Right + (a*(tao / (dy*dy)) - tao*Vcast / dy)*T_Left;
					}
					else if (j == 0 && i != 0 && i != nx - 1)//左边
					{
						T_New[i][j] = T_Cast;
					}
					else if (j == ny - 1 && i != 0 && i != nx - 1)//右边
					{
						T_Up = T_New[i + 1][j];
						T_Down = T_New[i - 1][j];
						T_Right = T_New[i][j - 1];
						T_Left = T_New[i][j - 1];
						T_Last[i][j] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + (1 - 2 * a*(tao / (dx*dx)) - 2 * a*(tao / (dy*dy)) + tao*Vcast / dy)*T_Last[i][j] +
							a*(tao / (dy*dy))*T_Right + (a*(tao / (dy*dy)) - tao*Vcast / dy)*T_Left;
					}
					else if (i == 0 && j == 0)//左下角
					{
						T_Last[i][j] = T_Cast;
					}
					else if (i == 0 && j == ny - 1)//右下角
					{
						h = Boundary_Condition(j, ny, Ly, ccml, H_Init);
						T_Up = T_New[i + 1][j];
						T_Down = T_New[i + 1][j] - 2 * dx*h*(T_New[i][j] - Tw) / lamda;
						T_Right = T_New[i + 1][j];
						T_Left = T_New[i][j - 1];
						T_Last[i][j] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + (1 - 2 * a*(tao / (dx*dx)) - 2 * a*(tao / (dy*dy)) + tao*Vcast / dy)*T_Last[i][j] +
							a*(tao / (dy*dy))*T_Right + (a*(tao / (dy*dy)) - tao*Vcast / dy)*T_Left;
					}
					else if (i == nx - 1 && j == 0)//左上角
					{
						T_Last[i][j] = T_Cast;
					}
					else if (i == nx - 1 && j == ny - 1)//右上角
					{
						T_Up = T_New[i - 1][j];
						T_Down = T_New[i - 1][j];
						T_Right = T_New[i - 1][j];
						T_Left = T_New[i][j - 1];
						T_Last[i][j] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + (1 - 2 * a*(tao / (dx*dx)) - 2 * a*(tao / (dy*dy)) + tao*Vcast / dy)*T_Last[i][j] +
							a*(tao / (dy*dy))*T_Right + (a*(tao / (dy*dy)) - tao*Vcast / dy)*T_Left;
					}
					else//内部
					{
						T_Up = T_New[i + 1][j];
						T_Down = T_New[i - 1][j];
						T_Right = T_New[i][j + 1];
						T_Left = T_New[i][j - 1];
						T_Last[i][j] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + (1 - 2 * a*(tao / (dx*dx)) - 2 * a*(tao / (dy*dy)) + tao*Vcast / dy)*T_Last[i][j] +
							a*(tao / (dy*dy))*T_Right + (a*(tao / (dy*dy)) - tao*Vcast / dy)*T_Left;
					}
				}
			}
		}
		
		if(k%out_num==0)
		{
			for(int i=0;i<nx;i++)			
				for(int j=0;j<ny;j++)
					T[i][j][count]=T_Last[i][j];			
			printf("Current time step=%d ",k);
			printf("Simulation time=%f ",k*tao);
			//Monitor_MeanTemperature(nx, ny, count, Y_Label);
			printf("\n");
			count++;
		}
		/*for(int i=0;i<nx;i++)		
			for(int j=0;j<ny;j++)	
				T_Last[i][j]=T_New[i][j];*/	
		flag = !flag;
	}

	end=clock();//结束计时
	printf("\n");
	printf("Assuming time =%d(ms)\n",(end-start));

	OutputTemperature(nx,ny,tnpts,out_num);

	return 0;
}
void Physicial_Parameters(float T)//固液参数
{
	float Ts = 1462.0, Tl = 1518.0, lamdas = 30, lamdal = 50, phos = 7000, phol = 7500, ce = 540.0, L = 265600.0, fs = 0.0;
	if (T<Ts)
	{
		fs = 0;
		pho = phos;
		lamda = lamdas;
		Ce = ce;
	}

	if (T >= Ts&&T <= Tl)
	{
		fs = (T - Ts) / (Tl - Ts);
		pho = fs*phos + (1 - fs)*phol;
		lamda = fs*lamdas + (1 - fs)*lamdal;
		Ce = ce + L / (Tl - Ts);
	}

	if (T>Tl)
	{
		fs = 1;
		pho = phol;
		lamda = lamdal;
		Ce = ce;
	}
}
float Boundary_Condition(int j,int ny,float Ly,float *ccml_zone,float *H_Init)
{
	float YLabel, h;
	YLabel = (j*Ly) / float(ny - 1);

	for (int i = 0; i < Section; i++)
	{
		if (YLabel >= *(ccml_zone + i) && YLabel <= *(ccml_zone + i + 1))
		{
			h = *(H_Init + i);
		}
	}
	return h;
}
void Monitor_MeanTemperature(int nx, int ny, int count, int *Y_Label)
{
	float temp_surface[Section] = { 0.0 }, temp_central[Section] = { 0.0 };
	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < Section; i++)
		{
			if (j >= *(Y_Label + i) && j <= *(Y_Label + i + 1))
			{
				temp_surface[i] = temp_surface[i] + T[0][j][count];
				temp_central[i] = temp_central[i] + T[nx-1][j][count];
			}
			printf("\n Surface temperature \n");
			float Mean_Temperature_Surface[Section] = { 0.0 }, Mean_Temperature_Central[Section] = { 0.0 };
			for (int i = 0; i < Section; i++)
			{
				Mean_Temperature_Surface[i] = temp_surface[i] / (*(Y_Label + i + 1) - *(Y_Label + i));
				printf("zone %d =%f  ", i + 1, Mean_Temperature_Surface[i]);
			}
			printf("\n Central temperature \n");
			for (int i = 0; i < Section; i++)
			{
				Mean_Temperature_Central[i] = temp_central[i] / (*(Y_Label + i + 1) - *(Y_Label + i));
				printf("zone %d =%f  ", i + 1, Mean_Temperature_Central[i]);
			}

		}
	}
}
void OutputTemperature(int nx,int ny,int tnpts,int out_num)
{
	FILE *fp=NULL;
	int k;

	fp=fopen("E:\\data_zf\\CPU2D.txt","w");
	k=(tnpts-1)/out_num-1;
	printf("Current time step is %d",k);
	for(int j=0;j<ny;j++)
	{
		for(int i=0;i<nx;i++)
		{
			fprintf(fp,"%f",T[i][j][k]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}