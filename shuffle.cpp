/*************************************
* THIS PROGRAM IS CALCULATING WITH A DIFFERENT "INITIAL CONDIT*
* ION". THE DIMERS ARE SPACED EVENLY WITH THE LOOPS EXTRACTED.* 
*                                                             *
* This file is to calculate the position of m dimer motors on *
* a lattice of size L.                                        *
* The position of each of the heads is stored in an array and *
* the calculation is done using the gillespie method.         *
*                                                             *
* The steps are:                                              *
* 1: The motors are displaced randomly along the lattice.     *
*                                                             *
* 2: For each head the rates for forward and backward motion  *
* are determined;                                             *
* left haeds: r_plus to left and r_minus to right (if empty)  *
* right heads: r_plus to right and r_minus to left (if empty) *
*                                                             *
* 3: The gillespie step: choose random number R and calculate *
* the time step; Delta_t = -1/r ln(R)                         *
*                                                             *
* 4: Process determination step:                              *
* Find non-zero r_plus and r_minus and index them.            *
* Choose another random variable p;                           *
* If p =< Sum(r_minus)/r ->choose another rondom number<index *
* Adjust the particle positions accordingly;                  *
* Similar steps for p > Sum(r_minus)/r                        *
*                                                             *
* 5: Take time t to t+Delta_t and repeat until t >= t_total   *         
***************************************************************
*                                                             * 
*/

#include <math.h>
#include <stdlib.h>
#include <fstream.h>
#include <time.h>
#include <iostream.h>
#include <unistd.h>


using namespace std;


//list of consts for ranq file from numerical recipes
//I'm taking the overflow at 2^31
#define IA 7141
#define IM 259200
#define IC 54773

// This structure is to be able to write the seed in a binary 
//file and be able to access it.




// model constants 
const int num =5;     //number of dimer motors
const int m = 10;      //number of motor heads; it should be 2*m
const int L = 200;   //lattice size
//const double t_total = 10000.0;
const double t_total = 2.0;
const float r_plus = 1.0;  //defining constants r_plus and r_minus
const float r_minus = 1.0;




/**********rate calculation functions********/ 

/*rate minus:look at the position of the particle behind for odd 
particles and particle ahead for even particles  */
//rates for the end particles should be calculated separately

void rates(int *lattice_ptr, int *position_ptr, int *rate_p_ptr, int *rate_m_ptr){
	
 	
	for (int i =0; i < m; i++){
		*(rate_m_ptr +i) = 0;
		*(rate_p_ptr +i) = 0;}
		
		
/* for plus particles r+ is going to lattice site ahead 
and r- is to go to lattice site behind; vice versa for minus particles         
RATES ARE PUT INTO LATTICE POINT CONTAINING PARTICLE AT BEGINNING OF THE MOVE.*/

	for (int i=1; i < (m-1); i++) {
		int x = *(position_ptr +i)-1;
		int latt_x = *(lattice_ptr +x);
		int latt_xp = *(lattice_ptr + x+1);
		int latt_xm = *(lattice_ptr + x-1);
		

		if ((latt_x > 0) && (latt_xp == 0)) 
			*(rate_p_ptr +i) = 1;
		else if ((latt_x > 0) && (latt_xm == 0))
			*(rate_m_ptr +i) = 1;
		else if ((latt_x < 0) && (latt_xp == 0))
			*(rate_m_ptr +i) = 1;
		else if ((latt_x < 0) && (latt_xm == 0))
			*(rate_p_ptr +i) = 1;
	}
	
	//end points
	
	if ((*position_ptr == 0) && (*(lattice_ptr +1) == 0))    //first lattice site
		*(rate_m_ptr) = 1;

	if (*(position_ptr + m-1) == L  && (*(lattice_ptr + L-2) == 0))  //last lattice site
		*(rate_m_ptr + L-1) = 1;
	

	return;
}








/**************gillespie step calculating the time interval for something to happen****/
 
double g_step(float r_total, unsigned long *jran_ptr){
	
	*jran_ptr = (*jran_ptr * IA + IC) % IM;
        double R = (double) *jran_ptr/(double) IM;
 	double z = log(R);
	double Delta_t = -1.0 /(r_total) * z;
	return(Delta_t);
}









/***********function determining if the event was r_plus*******/

int event_check(int minus, float r_total, unsigned long *jran_ptr){
	
        *jran_ptr = (*jran_ptr * IA + IC) % IM;
        float RR = (float) *jran_ptr/(float) IM;
	int check = 0;
			
	if ( RR <= ((minus*r_minus) /r_total))
		check = -1;
	else 
		check = +1;
	return(check);
}








/*******************procedure for shuffle*************/









/*******************procedure changing position for r_plus**********/

void event_forward(int *lattice_ptr,int *position_ptr, int *rate_p_ptr, unsigned long *jran_ptr){

	
	int test = 0;


	while (test ==0){
	
		*jran_ptr = (*jran_ptr * IA + IC) % IM;
       		int RRR =  m *(*jran_ptr) / IM;   //random motor head

		int y = *(position_ptr + RRR)-1;   //position of the chosen motor head
		int latt_y = *(lattice_ptr + y); //lattice site

		if ( (*(rate_p_ptr + RRR) != 0) && (latt_y < 0)){
			//cout <<RRR<<"	"<<"ratep"<< *(rate_p_ptr +RRR)<<"	";
			test++;
			*(lattice_ptr+y-1) = *(lattice_ptr+y);
			*(lattice_ptr+y) = 0;
			y--;
		}

		else if ((*(rate_p_ptr + RRR) != 0) && ((latt_y > 0) != 0)){
			//cout << "ratep" << *(rate_p_ptr+RRR) <<"	";
			*(lattice_ptr + y+1) = *(lattice_ptr+ y);
			*(lattice_ptr+y) = 0;
			y++;
			test++;
		}

		*(position_ptr+RRR) = y+1;
	}

	return;

}
	 









/***********************procedure changin position for r_minus**********/

void event_backward(int *lattice_ptr, int *position_ptr, int *rate_m_ptr, unsigned long *jran_ptr){


	int test =0;

	while (test == 0){
		
		*jran_ptr = (*jran_ptr * IA + IC) % IM;
                int RRR =  m *(*jran_ptr) / IM;

		int y = *(position_ptr + RRR)-1;
		int latt_y = *(lattice_ptr + y);

		if ((*(rate_m_ptr + RRR) != 0) && (latt_y < 0)){
			//cout << "ratem" << *(rate_m_ptr+RRR) <<"	";
			*(lattice_ptr +y+1)= *(lattice_ptr +y);
			*(lattice_ptr +y) = 0;
			y++;
			test++;
		}

		else if ((*(rate_m_ptr+ RRR) !=0 ) && (latt_y >0)){
			//cout <<*(rate_m_ptr+RRR) <<"	";
			*(lattice_ptr +y-1) = *(lattice_ptr + y);
			*(lattice_ptr + y) = 0;
			y--;
			test++;
		}

		*(position_ptr+RRR) = y+1;

	}

	return;


}











/***********************starting main*******************************/
/********************************************************************/




int main()
{ 
	int position[m];     //list of motor head positions
	int lattice[L];      //list of laatice sites
	int *position_ptr = &position[0];  /*pointer to the first poition */ 
	int *lattice_ptr = &lattice[0]; 
	int loop[num]; 
	
	int i,j,count;

	

//initiating

	for (i = 0; i < m; i++){
		position[i] = -1;}
	for (i = 0; i < L; i++){
		lattice[i] = 0;}
        
	
	//      opening the files 

	ofstream fout_position;  //declaring the file object
    fout_position.open("pm.txt", ios::app);                
	ofstream fout_loop;
	fout_loop.open("l.txt", ios::app);

	ifstream inData;
	inData.open("input.txt");
	inData >> i;
	

/**************seed initiation*********************/
	
	unsigned long ranj;
	unsigned long *jran_ptr;
	ranj = i;
    jran_ptr =&ranj;




/******************Initiation*************/    

	double t = 0.0;           //total time elapsed
//	fout_position << t <<"	";    //write the time in the file
//fout_loop << t <<"	";



//This is the normal initiation
	
/*	int counter = 0;

	while (counter < num){

	    *jran_ptr = (*jran_ptr * IA + IC) % IM;
	    int x =  L *(*jran_ptr) / IM;
	    if ( (lattice[x-1] == 0) && (lattice[x] == 0) ){
		position[2*counter] = x;
		position[(2*counter + 1)] = x+1;
		lattice[x-1] = -(counter+1);
		lattice[x] = counter+1;
		counter++;}
	}	
*/

// DIFFFERENT INITIAL CONDITIONS


	for (int count = 0; count < num; i++){
		int x = (int) L/num * count +1; 
		position[(2*count)] = x;
		lattice[x-1]= -(count+1);
		int y = (int) L/num *(count+1)-1;
		position[2*count+1]= y;
		lattice[y-1] = count+1;
	}

	
	for (int count = 0; count < m; count++){
		int proxy = position[count];
		fout_position << position[count]<<"	"<< lattice[proxy - 1]<<"	";}
	fout_position << "\n";
          



//end of initiation cycle

/*******************************************************/

	int rate_p[m];
	int *rate_p_ptr = &rate_p[0];

	int rate_m[m];
	int *rate_m_ptr = &rate_m[0];	
	//int mycounter = 0;



//recursion cycle to go through t_total in delta_t steps; each delta_t comes from the gillespie cycle

        //while (t <= t_total){
			
		int plus = 0;  //number of particles with non-zeros r_plus and r_minus
		int minus = 0;
			

		for (i =0; i < m; i++){
			*(rate_m_ptr+i) = 0;
			*(rate_p_ptr+i) = 0;}

		rates( lattice_ptr, position_ptr, rate_p_ptr, rate_m_ptr);
			
		for (int n = 0;  n < m; n++){
			if (*(rate_m_ptr + n) != 0)
				minus++;
			if (*(rate_p_ptr+n) != 0)
				plus++;}
		float r_total = plus * r_plus + minus* r_minus; //total rate




		double delta_t = g_step(r_total, jran_ptr);   //recalling the gillespie step 
		t = t + delta_t;

		int s = event_check(minus, r_total, jran_ptr);
		if (s == -1)
			event_backward(lattice_ptr, position_ptr,rate_m_ptr, jran_ptr);
		if (s == +1)
			event_forward(lattice_ptr, position_ptr, rate_p_ptr, jran_ptr);

		

//		fout_position << t <<"	"<< s <<"\n";
		//fout_loop << t << "	";
                        

			
		for (int count = 0; count < num; count++) {
			loop[count] = position[(2*count+1)] - position[(2*count)] + 1;
	//	fout_position << position[2*count] <<"	" << position[2*count +1] << "	";
			//fout_loop << loop[count] <<"\n";
		}
		
	for (int myi = 0; myi < m; myi++)
		{int myy= position[myi]-1;
		fout_position << position[myi]<<"	"<<lattice[myy]<<"	";}

	fout_position << "\n";
		//fout_loop <<"\n";		
//	} 
       

 
	for (int count = 0;  count < num; count++){ 
		fout_loop << loop[count] << "\n";}
	//fout_loop <<"\n";


	ofstream outData;
	outData.open("input.txt");

	outData << ranj <<"\n";
	outData.close();
	inData.close();


	fout_position.close();
	fout_loop.close();
	
	return(0);
}               


















                                      
