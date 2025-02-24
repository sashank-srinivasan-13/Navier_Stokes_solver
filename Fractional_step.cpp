#include <iostream>
#include <fstream>
#include <cmath>

int Ny, Nx;
double **P, **P_old, **Sm, **u, **v, **u_hat, **v_hat, **u_old, **v_old;
double dx,dy,dx2,dy2;
int QUICK=1; //0 for Upwind, 1 for QUICK, -1 For central fluxes

double QUICKuxp(int i , int j)
{
	double uret=0;

	if(j==0)
		uret=(u_old[i][j]+u_old[i][j+1])/2.0;
	else if(QUICK==1)
	{
		uret=6.0*u_old[i][j]+3.0*u_old[i][j+1]-1.0*u_old[i][j-1];
		uret=uret/8.0;
	}
	else
		uret=u_old[i][j];


	return uret;
}

double QUICKuxn(int i , int j)
{
	double uret=0;

	if(j==Nx-2)
		uret=(u_old[i][j]+u_old[i][j-1])/2.0;
	else if(QUICK==1)
	{
		uret=6.0*u_old[i][j+1]+3.0*u_old[i][j]-1.0*u_old[i][j+2];
		uret=uret/8.0;
	}
	else
		uret=u_old[i][j+1];


	return uret;
}

double QUICKuyp(int i, int j)
{
	double uret=0;

	if(i==0)
		uret=(u_old[i][j]+u_old[i+1][j])/2.0;
	else if(QUICK==1)
	{
		uret=6.0*u_old[i][j]+3.0*u_old[i+1][j]-1.0*u_old[i-1][j];
		uret=uret/8.0;
	}
	else
		uret=u_old[i][j];


	return uret;
}

double QUICKuyn(int i, int j)
{
	double uret=0;

	if(i==Ny-2)
		uret=(u_old[i][j]+u_old[i-1][j])/2.0;
	else if(QUICK==1)
	{
		uret=6.0*u_old[i+1][j]+3.0*u_old[i][j]-1.0*u_old[i+2][j];
		uret=uret/8.0;
	}
	else
		uret=u_old[i+1][j];


	return uret;
}


double QUICKvxp(int i , int j)
{
	double vret=0;

	if(j==0)
		vret=(v_old[i][j]+v_old[i][j+1])/2.0;
	else if(QUICK==1)
	{
		vret=6.0*v_old[i][j]+3.0*v_old[i][j+1]-1.0*v_old[i][j-1];
		vret=vret/8.0;
	}
	else
		vret=v_old[i][j];


	return vret;
}

double QUICKvxn(int i , int j)
{
	double vret=0;

	if(j==Nx-2)
		vret=(v_old[i][j]+v_old[i][j-1])/2.0;
	else if(QUICK==1)
	{
		vret=6.0*v_old[i][j+1]+3.0*v_old[i][j]-1.0*v_old[i][j+2];
		vret=vret/8.0;
	}
	else
		vret=v_old[i][j+1];


	return vret;
}


double QUICKvyp(int i, int j)
{
	double vret=0;

	if(i==0)
		vret=(v_old[i][j]+v_old[i+1][j])/2.0;
	else if(QUICK==1)
	{
		vret=6.0*v_old[i][j]+3.0*v_old[i+1][j]-1.0*v_old[i-1][j];
		vret=vret/8.0;
	}
	else
		vret=v_old[i][j];

	return vret;
}

double QUICKvyn(int i, int j)
{
	double vret=0;

	if(i==Ny-2)
		vret=(v_old[i][j]+v_old[i-1][j])/2.0;
	else if(QUICK==1)
	{
		vret=6.0*v_old[i+1][j]+3.0*v_old[i][j]-1.0*v_old[i+2][j];
		vret=vret/8.0;
	}
	else
		vret=v_old[i+1][j];


	return vret;
}

double calc_Cu(int i, int j)
{
	double Cu=0;

	double a,b,c,d;

	if(u_old[i][j]>0)
	{
		a=QUICKuxp(i,j);
		b=QUICKuxp(i,j-1);

	}
	else
	{

		a=QUICKuxn(i,j);
		b=QUICKuxn(i,j-1);

	}


	if(v_old[i][j]>0)
	{
		c=QUICKuyp(i,j);
		d=QUICKuyp(i-1,j);
	}
	else
	{

		c=QUICKuyn(i,j);
		d=QUICKuyn(i-1,j);

	}


	if(QUICK==-1)
	{
		a=(u_old[i][j]+u_old[i][j+1])/2.0;
		b=(u_old[i][j]+u_old[i][j-1])/2.0;
		c=(u_old[i][j]+u_old[i+1][j])/2.0;
		d=(u_old[i][j]+u_old[i-1][j])/2.0;
	}

	Cu=u_old[i][j]*(a-b)/dx+v_old[i][j]*(c-d)/dy;


	return Cu;
}

double calc_Cv(int i, int j)
{
	double Cv=0;

	double a,b,c,d;


	if(u_old[i][j]>0)
	{
		a=QUICKvxp(i,j);
		b=QUICKvxp(i,j-1);
	}
	else
	{

		a=QUICKvxn(i,j);
		b=QUICKvxn(i,j-1);

	}

	if(v_old[i][j]>0)
	{
		c=QUICKvyp(i,j);
		d=QUICKvyp(i-1,j);
	}
	else
	{

		c=QUICKvyn(i,j);
		d=QUICKvyn(i-1,j);

	}


	if(QUICK==-1)
	{
		a=(v_old[i][j]+v_old[i][j+1])/2.0;
		b=(v_old[i][j]+v_old[i][j-1])/2.0;
		c=(v_old[i][j]+v_old[i+1][j])/2.0;
		d=(v_old[i][j]+v_old[i-1][j])/2.0;
	}

	Cv=u_old[i][j]*(a-b)/dx+v_old[i][j]*(c-d)/dy;

	return Cv;
}

double calc_Du(int i, int j)
{
	double Du=0;

	Du=(u_old[i+1][j]+u_old[i-1][j]-2.0*u_old[i][j])/dy2+(u_old[i][j+1]+u_old[i][j-1]-2.0*u_old[i][j])/dx2;

	return Du;
}

double calc_Dv(int i, int j)
{
	double Dv=0;

	Dv=(v_old[i+1][j]+v_old[i-1][j]-2.0*v_old[i][j])/dy2+(v_old[i][j+1]+v_old[i][j-1]-2.0*v_old[i][j])/dx2;

	return Dv;
}



int main(void)
{

	std::cout<<"Program started"<<std::endl;

	int i,j;
	int t=0;
	int iter_P;
	double omega;
	double Lx,Ly;
	double x,y;
	double Res_P, Res_global, Res_threshold;

	int nsteps;
	double dt,dt_limit;
	double time,t_final;

	double rho, mu, Re, V_boundary;

	FILE *fpt ;
	fpt = fopen ("Time_series_solution.dat" , "w+" );
	fprintf(fpt,"iteration time time_step_limit time_step_used Residual Pressure_solver_iterations \n");



	Ny=42; //40 grid points + 2 ghost nodes
	Nx=42; //40 grid points + 2 ghost nodes
	omega=1.4;
	Lx=1.0;
	Ly=1.0;
	dy=Ly/(Ny-2);
	dx=Lx/(Nx-2);
	dx2=dx*dx;
	dy2=dy*dy;

	t_final=10.0;
	time=0;
	nsteps=int(t_final/dt);
	Res_global=0;
	Res_threshold=5*pow(10,-5);

	Re=50;
	rho=1.0;
	V_boundary=1.0;
	mu=rho*V_boundary*Lx/Re;


	P=new double* [Ny];
	P_old=new double* [Ny];
	Sm=new double* [Ny];

	u=new double* [Ny];
	u_hat=new double* [Ny];
	u_old=new double* [Ny];
	v=new double* [Ny];
	v_hat=new double* [Ny];
	v_old=new double* [Ny];

	for(i=0;i<Ny;i++)
	{
		P[i]=new double [Nx];
		P_old[i]=new double [Nx];
		Sm[i]=new double [Nx];


		u[i]=new double [Nx];
		u_hat[i]=new double [Nx];
		u_old[i]=new double [Nx];
		v[i]=new double [Nx];
		v_hat[i]=new double [Nx];
		v_old[i]=new double [Nx];
	}


	//Initialization
	{

		for(i=0;i<Ny;i++)
			for(j=0;j<Nx;j++)
			{
				P[i][j]=pow(10,-10);
				P_old[i][j]=pow(10,-10);

				u[i][j]=pow(10,-10);
				v[i][j]=pow(10,-10);

				u_hat[i][j]=pow(10,-10);
				v_hat[i][j]=pow(10,-10);

				u_old[i][j]=pow(10,-10);
				v_old[i][j]=pow(10,-10);
			}

	}

	std::cout<<"Starting time loop"<<std::endl;


	//for(t=0;t<nsteps;t++)
	t=0;
	while(time<t_final)
	{

		t++;

		//Update values
		{
			for(i=1;i<Ny-1;i++)
				for(j=1;j<Nx-1;j++)
				{
					u_old[i][j]=u[i][j];
					v_old[i][j]=v[i][j];
					P_old[i][j]=P[i][j];
				}
		}

		//Calculate time step
		{
			dt_limit=0.01;
			double dummy;

			for(i=1;i<Ny-1;i++)
				for(j=1;j<Nx-1;j++)
				{
					dummy=fabs(u_old[i][j])/dx+fabs(v_old[i][j])/dy+2*mu/rho*(1/dx2+1/dy2);
					dummy=1.0/dummy;
					if(dummy<dt_limit)
						dt_limit=dummy;
				}

			dt=5*pow(10,-3);

		}

		std::cout<<t<<" "<<time<<" "<<dt_limit<<" "<<Res_global<<std::endl;


		//Velocity boundary conditions
		{
			j=0;
			for(i=1;i<Ny-1;i++)
			{
				u_old[i][j]=-u_old[i][j+1];
				v_old[i][j]=-v_old[i][j+1];
			}

			j=Nx-1;
			for(i=1;i<Ny-1;i++)
			{
				u_old[i][j]=-u_old[i][j-1];
				v_old[i][j]=-v_old[i][j-1];
			}

			i=0;
			for(j=1;j<Nx-1;j++)
			{
				u_old[i][j]=-u_old[i+1][j];
				v_old[i][j]=-v_old[i+1][j];
			}

			i=Ny-1;
			for(j=1;j<Nx-1;j++)
			{
				u_old[i][j]=2*V_boundary-u_old[i-1][j];
				v_old[i][j]=-v_old[i-1][j];
			}
		}

		//Calculate fractional velocities
		{
			double Cu,Cv,Du,Dv,Hu,Hv;


			for(i=1;i<Ny-1;i++)
				for(j=1;j<Nx-1;j++)
				{

					Cu=calc_Cu(i,j);
					Cv=calc_Cv(i,j);
					Du=calc_Du(i,j);
					Dv=calc_Dv(i,j);

					Hu=-Cu+mu/rho*Du;
					Hv=-Cv+mu/rho*Dv;

					u_hat[i][j]=u_old[i][j]+dt*(Hu);
					v_hat[i][j]=v_old[i][j]+dt*(Hv);

				}




		}

		//HatVelocity boundary conditions
		{
			j=0;
			for(i=1;i<Ny-1;i++)
			{
				u_hat[i][j]=-u_hat[i][j+1];
				v_hat[i][j]=-v_hat[i][j+1];
			}

			j=Nx-1;
			for(i=1;i<Ny-1;i++)
			{
				u_hat[i][j]=-u_hat[i][j-1];
				v_hat[i][j]=-v_hat[i][j-1];
			}

			i=0;
			for(j=1;j<Nx-1;j++)
			{
				u_hat[i][j]=-u_hat[i+1][j];
				v_hat[i][j]=-v_hat[i+1][j];
			}

			i=Ny-1;
			for(j=1;j<Nx-1;j++)
			{
				u_hat[i][j]=2*V_boundary-u_hat[i-1][j];
				v_hat[i][j]=-v_hat[i-1][j];
			}
		}

		//Calculate pressure source term from mass residual
		{
			double a,b,c,d;

			for(i=1;i<Ny-1;i++)
				for(j=1;j<Nx-1;j++)
				{

					a=(u_hat[i][j+1]+u_hat[i][j])/2.0;
					b=(u_hat[i][j-1]+u_hat[i][j])/2.0;
					c=(v_hat[i+1][j]+v_hat[i][j])/2.0;
					d=(v_hat[i-1][j]+v_hat[i][j])/2.0;



					Sm[i][j]=rho/dt*((a-b)/dx+(c-d)/dy);

				}

		}


		//Pressure calculation procedure
		{

			int flag=1;
			iter_P=0;


			while(flag)
			{


				//Pressure boundary conditions
				{
					j=0;
					for(i=1;i<Ny-1;i++)
						P[i][j]=P_old[i][j+1];

					j=Nx-1;
					for(i=1;i<Ny-1;i++)
						P[i][j]=P_old[i][j-1];

					i=0;
					for(j=1;j<Nx-1;j++)
						P[i][j]=P_old[i+1][j];

					i=Ny-1;
					for(j=1;j<Nx-1;j++)
						P[i][j]=P_old[i-1][j];
				}


				//Solve for pressure using SOR
				{



					for(i=1;i<Ny-1;i++)
						for(j=1;j<Nx-1;j++)
						{

							P[i][j]=
									+1.0/dx2*(P[i][j+1]+P[i][j-1])
									+1.0/dy2*(P[i+1][j]+P[i-1][j])
									-Sm[i][j];

							P[i][j]=P[i][j]/(2.0/dx2+2.0/dy2);

							P[i][j]=P_old[i][j]+omega*(P[i][j]-P_old[i][j]);


						}

				}

				//cacluclate Residual
				{
					double Rij=0.0;

					Res_P=0.0;

					for(i=1;i<Ny-1;i++)
						for(j=1;j<Nx-1;j++)
						{

							Rij=
									+(P[i][j+1]+P[i][j-1]-2.0*P[i][j])/dx2
									+(P[i+1][j]+P[i-1][j]-2.0*P[i][j])/dy2
									-Sm[i][j];

							Res_P=Res_P+Rij*Rij;
						}


					Res_P=Res_P/((Ny-2)*(Nx-2));
					Res_P=pow(Res_P,0.5);
				}

				if(Res_P<Res_threshold)
					flag=0;

				iter_P++;

				for(i=1;i<Ny-1;i++)
					for(j=1;j<Nx-1;j++)
						P_old[i][j]=P[i][j];

			}

		}


		//Update velocities at next time step
		{
			for(i=1;i<Ny-1;i++)
				for(j=1;j<Nx-1;j++)
				{

					u[i][j]=u_hat[i][j]-dt/rho*(P[i][j+1]-P[i][j-1])/(2.0*dx);
					v[i][j]=v_hat[i][j]-dt/rho*(P[i+1][j]-P[i-1][j])/(2.0*dy);

				}

		}

		time=time+dt;

		Res_global=0;
		for(i=1;i<Ny-1;i++)
			for(j=1;j<Nx-1;j++)
			{
				Res_global=Res_global+pow(u[i][j]-u_old[i][j],2.0)+pow(v[i][j]-v_old[i][j],2.0);
			}
		Res_global=pow(Res_global,0.5)/((Ny-2)*(Nx-2));

		fprintf(fpt ,"%d %e %e %e %e %d \n", t, time, dt_limit, dt , Res_global, iter_P);



		if(t*dt==0.5||t*dt==1.0||t*dt==2.0||t*dt==10.0)
		{
			FILE *fp ;
			char filename[250];


			sprintf(filename,"Solution_%d.dat",t);
			fp = fopen (filename , "w+" );
			fprintf(fp,"x y u v P \n");

			for(i=1;i<Ny-1;i++)
				for(j=1;j<Nx-1;j++)
				{
					x=(j-1)*dx+dx/2;
					y=(i-1)*dy+dy/2;
					fprintf(fp ,"%f %f %f %f %f \n", x, y, u[i][j], v[i][j], P[i][j]);
				}
			fclose(fp);

			sprintf(filename,"u_Centerline_%d.dat",t);
			fp = fopen (filename , "w+" );
			fprintf(fp,"y u \n");

			for(i=1;i<Ny-1;i++)
			{
				y=(i-1)*dy+dy/2;
				fprintf(fp ,"%f %f \n",y, (u[i][Nx/2-2]+u[i][Nx/2-1])/2);
			}
			fclose(fp);

			sprintf(filename,"v_Centerline_%d.dat",t);
			fp = fopen (filename , "w+" );
			fprintf(fp,"x v \n");

			for(j=1;j<Ny-1;j++)
			{
				x=(j-1)*dx+dx/2;
				fprintf(fp ,"%f %f \n",x, (v[Ny/2-2][j]+v[Ny/2-1][j])/2);
			}
			fclose(fp);
		}

	}

	fclose(fpt);


	FILE *fp ;
	fp = fopen ("Final_Solution.dat" , "w+" );
	fprintf(fp, "Lx = %f Ly = %f Reynolds number = %f Density = %f Viscosity = %f Final time = %f \n", Lx, Ly, Re, rho, mu, time);

	fprintf(fp,"x y u v P \n");
	for(i=1;i<Ny-1;i++)
		for(j=1;j<Nx-1;j++)
		{
			x=(j-1)*dx+dx/2;
			y=(i-1)*dy+dy/2;
			fprintf(fp ,"%f %f %f %f %f \n", x, y, u[i][j], v[i][j], P[i][j]);
		}

	fclose(fp);

	for(i=0;i<Ny;i++)
	{
		delete [] P[i];
		delete [] P_old[i];
		delete [] Sm[i];

		delete [] u_hat[i];
		delete [] u_old[i];
		delete [] u[i];

		delete [] v_hat[i];
		delete [] v_old[i];
		delete [] v[i];
	}

	delete [] P;
	delete [] P_old;
	delete [] Sm;

	delete [] u_hat;
	delete [] u_old;
	delete [] u;

	delete [] v_hat;
	delete [] v_old;
	delete [] v;

	return 0;
}


