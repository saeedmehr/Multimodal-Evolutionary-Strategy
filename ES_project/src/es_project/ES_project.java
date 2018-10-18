
package es_project;

import java.util.Scanner;
import java.util.Random;
public class ES_project {

    public static  int mu;   //valedin
    public static  int lamda; // farzandan
    public static  int dimension;
    public static  int alphasize;
    static int [][] interval;
    static double stepsize=2;
    static double Tau1 = 0.15;
    static double Tau2 = 0.15;
    static double betha = 0.0873;
    //static double  LowerLimitStepSize  = 0.0000005;
static double  LowerLimitStepSize  = 0.0001;


    static Individual [] tpop;        //temprory population
    public static void main(String[] args) {
        Random r = new Random ();
        Individual best=new Individual();
        best.fitness=Double.MAX_VALUE;
        Scanner sc= new Scanner(System.in);
        System.out.println("please insert your population size:\n");
        mu=sc.nextInt();
        lamda=7*mu;  //halate behine.jozve 46
        tpop=new Individual [lamda];
        System.out.println("please insert your dimension size:\n");
        dimension=sc.nextInt();
        alphasize=dimension*(dimension-1)/2;  //kolitarin halate momken.43
        interval= new int[dimension][2];
        for (int i = 0; i < dimension; i++) {
            System.out.println("enter the "+(i+1)+"'th lower bound:\n");
            interval[i][0]=sc.nextInt();
            System.out.println("enter the "+(i+1)+"'th upper bound:\n");
            interval[i][1]=sc.nextInt();

        }

        Individual [] i = new Individual[mu];    //mu ta individual vase jamiate avalie
        for (int j = 0; j < mu; j++) {
            i[j]=new Individual();
            for (int k = 0; k < dimension; k++) {
                i[j].input[k]=interval[k][0]+r.nextDouble()*(interval[k][1]-interval[k][0]); //random beine kamtarin ta bishtarin
                i[j].step[k]=stepsize;
            }       // ta inja population haye avalie ro dar ovordim


            for (int a = 0; a < alphasize; a++) {
                i[j].angle[a]=(r.nextDouble()-0.5)*2*Math.PI;  // alpha =>[-pi , pi].46


            }
            fitFunc(i[j]);
            if(i[j].fitness<best.fitness)
                best=i[j];

        }
        /** #################################################################################### **/
        Individual parent1;
        Individual parent2;
        Individual child;
        int naslindex=0;
        for (int j = 0; j <1000; j++) {   //naasl

            //recombination uniform 
            for (int k = 0; k < lamda; k++) {
                child=new Individual();

                for (int l = 0; l < dimension; l++) {// x hara crossoverh kardim faghat ta inja 
                    parent1 = i[r.nextInt(mu)];
                    parent2 = i[r.nextInt(mu)];
                    if(r.nextBoolean())
                        child.input[l]=parent1.input[l];
                    else
                        child.input[l]=parent2.input[l];

                    child.step[l]=(parent1.step[l]+parent2.step[l])/2;   //cross overe e sigma ha.46(baztarkibie global)

                }
                for (int l = 0; l < alphasize; l++) {
                    parent1 = i[r.nextInt(mu)];
                    parent2 = i[r.nextInt(mu)];
                    child.angle[l]=(parent1.angle[l]+parent2.angle[l])/2;
                }

                //Mutate Sigmas
                double tmpNoise=r.nextGaussian();
                for (int q = 0; q <dimension; q++) {
                    child.step[q] *= Math.exp(Tau1* tmpNoise + Tau2*r.nextGaussian()); //45. ba tavajoh be jozve
                    if (child.step[q] < LowerLimitStepSize)
                        child.step[q] =LowerLimitStepSize;
                }
                //Mutate Alphas.45
                for (int q = 0; q < alphasize; q++) {
                    child.angle[q] += betha*r.nextGaussian();
                    if ( child.angle[q]< -Math.PI)
                        child.angle[q] = -Math.PI;
                    if (child.angle[q] > Math.PI)
                        child.angle[q] = Math.PI;
                }
                double[]  xCopy = new double[dimension]; //44safhe akhar
                //Generate mutation vector in unitspace modified by sigmas
                for (int q = 0; q < dimension; q++) {
                    xCopy[q] = child.step[q]*r.nextGaussian();     //44safhe akhar.jozve
                }

                //turn mutationvector with alphas
                for (int q = 0; q < dimension-1; q++) { //az n-1
                    for (int z = q+1; z < dimension; z++) { // ta n
                        double alpha=getAlpha(q,z,dimension,child);
                        double xX=java.lang.Math.cos(alpha)*xCopy[q]-java.lang.Math.sin(alpha)*xCopy[z]; //batavajoh be 48
                        double xY=java.lang.Math.sin(alpha)*xCopy[q]+java.lang.Math.cos(alpha)*xCopy[z]; //batavajoh be 48
                        xCopy[q]=xX;
                        xCopy[z]=xY;
                    }
                }

                //modify genotype
                for (int q = 0; q < dimension; q++) {
                    child.input[q] += ((interval[q][1] -interval[q][0])/2)*xCopy[q];
                    if (interval[q][0] > child.input[q])
                        child.input[q] = interval[q][0];
                    if (interval[q][1] < child.input[q])
                        child.input[q] = interval[q][1];
                }
                fitFunc(child);
                tpop[k]=child;



            }
            Individual temp;
            for(int q=0; q < lamda; q++){
                for(int m=1; m < (lamda-q); m++){

                    if(tpop[m-1].fitness > tpop[m].fitness){
                        //swap the elements!
                        temp = tpop[m-1];
                        tpop[m-1] = tpop[m];
                        tpop[m] = temp;
                    }

                }
            }
            for (int k = 0; k < 10; k++) {
                i[k]=tpop[k];
            }

            if(best.fitness>i[0].fitness) {
                best = i[0];

            naslindex=j;
            }



        }
        System.out.println("The best Individual fitness was: "+best.fitness+"\n"
                + "Inputs : \n"+toStr(best.input)+"\ngeneration number:"+naslindex);





    }

    public static double getAlpha(int i, int j, int n,Individual child) {
        int sum=0;
        if (i<j) {
            for (int count=0; count<i; count++) {
                sum+=n-count-1;
            }
            sum+=j-i;
            sum--;
            return child.angle[sum];
        }else{
            return 0.0; //wrong access to alpha list
        }
    }

    public static String toStr(double[] input) {
        String res="";
        for (int i = 0; i < input.length; i++) {

            res+="input "+i+"'th"+input[i]+" ";
        }
        return res;
    }

    // ma az in class be tedade muta new mikonim bad 7muta ke to jozve hast tozish
    public static class Individual{
        double [] input = new double[dimension];         // x haye mas darvaghe
        double [] step = new double[dimension];           //hamon sigmahas
        double [] angle = new double[dimension*(dimension-1)/2];           //hamon alphst
        double  fitness;

    }
    public static void fitFunc(Individual p){ //person as individual


        //*************************************f2***************************************
   double sum1=0;
    double sum2=1;
    int dimension=p.input.length;
    for (int i = 0; i < dimension; i++)
    {
      sum1 = sum1 + p.input[i] * p.input[i];
      sum2*=Math.sin(p.input[i]) ;
    }
    p.fitness =-0.0001*Math.pow((Math.abs(sum2*Math.exp(Math.abs(100-Math.sqrt(sum1)/Math.PI)))+1),0.1);
    }}
//******************************************f3********************************************

//        double sum1=0;
//        double sum2=0;
//        int dimension=p.input.length;
//        for (int i = 0; i < dimension; i++)
//        {
//            if(i==0){
//                sum1=(p.input[0]-1)*(p.input[0]-1);
//            }else{
//                sum2+=(i+1)*(2*p.input[i]*p.input[i]-p.input[i-1])*(i+1)*(2*p.input[i]*p.input[i]-p.input[i-1]);
//            }}
//        p.fitness =sum1+sum2;
//    }}

//***************************************F1************************************************    
//     double sum1=0;
//    double sum2=0;
//        for (int i = 0; i < dimension; i++)
//{
//sum1 = sum1 + p.input[i] * p.input[i];              //gen hamon inpute mane
//sum2 = sum2 + (Math.cos(2 * Math.PI * p.input[i]));  //to java cos o exp kochike to c# bozorge
//}
//p.fitness = -20 * Math.exp(-0.2 * Math.sqrt(sum1 / dimension)) - Math.exp(sum2 / dimension) + 20 + Math.exp(1);
//    }}

//*************************tabeye digar******************************//
//        int dimension=p.input.length;
//    double sum= 418.9829 * dimension;
//    for (int j = 0; j < dimension; j++) {
//            if (Math.abs(p.input[j]) > 500)
//       p.fitness= Double.MAX_VALUE;             
//      sum -= p.input[j]*Math.sin(Math.sqrt(Math.abs(p.input[j])));
//        }
//        p.fitness=sum;
//        //*******************************************************************

// }

