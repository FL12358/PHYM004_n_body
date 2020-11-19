#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>

/*  version     =   1.5.2
    author      =   Frederick W Liardet
    date        =   2019.12.13


    nbody.c takes a data file in a tab separated format, this file contains the name of the body, mass, position (3D, in km) and
    velocity (3D, in m/s). Comments at the start of the file are designated with a line starting with '#'. An example file
    (OrbitDataESM.txt) has the layout of:
                # Name	Mass	x	y	z	Vx	Vy	Vz
                Earth	5.9726e24	1.5e8	0.0	0.0	0.0	2.929e4	0.0
                Sun	1.9885e30	0.0	0	0.0	0	0.0	0.0
                Moon	7.342e22	1.50363e8	0.0	0.0	0.0	3.0312e4	0.0

    Care has to be made to make sure that the file ends with a new line character '\n', on a windows machine this may not be
    added automatically and will cause errors.

    The command line arguments are as follows:
        -r  Runs the simulation with the 4th order Runge Kutta method, followed by the data file name
        -v  Runs the simulation with the Velocity Verlet method, follows by the data file name
        -s  Total time span of the simulation, units of seconds, double float format(Defaults to STD_TIME_LEN)
        -t  Time step between each iteration of the simulation, units of seconds, double float format (Defaults to STD_STEP)

    Example compile command:
        $ gcc -Wall nbody.c -lm -lgsl > nbodyProg

    Note: The -lgslcblas argument may also be needed to compile the program
    Example simulation command (with a time span of 100 days and a time step of 120 seconds):
        $ ./nbodyProg -v OrbitDataSolSys.txt -s 8640000 -t 120

    Command line output from command:
        $   Method:     Velocity Verlet
        $   Number of bodies detected:  9
        $       (Earth, Venus, Mars, Mercury, Jupiter, Saturn, Uranus, Neptune, Sun)
        $   Orbital Periods Retrieved (Days):
        $       Mercury:         88.76
        $   Total Energy (J):
        $       Initial:         -1.965E+35
        $       Final:           -1.965E+35
        $       Change:          1.496E+30
        $       Percent Error:   0.000761%
        $   Execution Time (s):      8.9220s

    This will also create a text file with name OUT_FILENAME, which has the position values for all the bodies in a 3 long tab
    separated format (this loses the information about the differentiation between the bodies, however it will plot fine). A
    command is also run using GNUPlot and the 'nBodyPlot.script' file to plot the paths of the bodies in the x-y plane, this
    plot is then saved to 'nBodyImage.png'.

*/
#define GNUPLOT_EXE    "gnuplot"
#define GNUPLOT_SCRIPT "nBodyPlot.script"
#define MAX_LINE_SIZE 1000
#define MAX_COMMAND_LEN 1000
#define DIMS 3                              //  Dimensions of system
#define BODY_REF_FRAME -1                   //  Body index of reference frame for plotting (-1 is based on static coordinates)
#define FILE_PRINT_NO 100                   //  Prints to file every FILE_PRINT_NO iteration

static const double BIG_G = 6.67430e-11;    //Standard value of G
static const double SOFT_L = 1e2; //  in m
static const double STD_STEP = 60;
static const double STD_TIME_LEN = 1.0e7;
static const char OUT_FILENAME[MAX_LINE_SIZE] = "OUTPUT.txt";

typedef struct{
    double pos[DIMS]; //[x,y,z...] positions m
    double vel[DIMS]; //[x,y,z...] velocities m/s
    double acc[DIMS]; //[x,y,z...] accelerations m/s^2
    double mass; //Mass in kg
    char name[MAX_LINE_SIZE]; //Name of body
} BODY;

typedef struct{ //Used to pass parameters into RK4 function
    BODY* bodies;
    int bodyNo;
    int bodyID;
} PARAMS;

void* xmalloc(size_t bytes) {
/*  functionally the same as malloc(), but throws error when memory allocation fails
*/
    void* retVal = malloc(bytes);
    if (retVal) {
        return retVal;
    }
    fprintf(stderr,"ERROR: Failed to allocate memory");
    exit(EXIT_FAILURE);
}

static double ReadFileBodyNo(char* filename){
/*  Reads in file and returns the number of bodies in the file
    Returns -1 if file does not exist
    Returns 0 if no bodies are in file
*/
    FILE *file;
    char line[MAX_LINE_SIZE];
    char temp;
    int bodyNo = 0;
    file = fopen(filename, "r");

    if(file == NULL){
        return -1;  // Error reading file
    }
    while (1){  //Gets the dimensions of the matrix, skips all rows starting with "#"
        fgets(line, MAX_LINE_SIZE, file);
        if (line[0] != '#'){
            break;
        }
    }
    while(temp!=EOF){
        temp = getc(file);
        if(temp == '\n'){
            bodyNo = bodyNo + 1;
        }
    }
    fclose(file);
    return bodyNo + 1;
}

static BODY* ReadFileData(char* filename, int bodyNo){
/*  Reads in file and assigns the data from the file to a 'bodyNo' long array of BODY structs
    Converts the distances in the data file into meters and sets the acceleration values to zero
    Returns the array of BODY structs
*/
    BODY* bodies = xmalloc(bodyNo*sizeof(BODY));
    char line[MAX_LINE_SIZE];
    FILE *file;
    file = fopen(filename, "r");
    int scanNo = 0;

    while (1){  //Gets the dimensions of the matrix, skips all rows starting with "#"
        fgets(line, MAX_LINE_SIZE, file);
        if (line[0] != '#'){
            break;
        }
    }
    for(int i=0;i<bodyNo;i++){
        scanNo = sscanf(line, "%s %lf %lf %lf %lf %lf %lf %lf",
            bodies[i].name, &bodies[i].mass,
            &bodies[i].pos[0], &bodies[i].pos[1], &bodies[i].pos[2],
            &bodies[i].vel[0], &bodies[i].vel[1], &bodies[i].vel[2]);
        if(scanNo != 8){    //  8 number of data pieces per body
            fprintf(stderr,"ERROR: Could not read in data from file\n");
            exit(EXIT_FAILURE);
        }
        fgets(line, MAX_LINE_SIZE, file);
    }
    for(int i=0;i<bodyNo;i++){
        for(int j=0;j<DIMS;j++){
        bodies[i].pos[j]*=1000;
        bodies[i].acc[j] = 0;
        }
    //printf("NAME: %s\n",bodies[i].name);
    }
    return bodies;
}

static void GnuPlot(void){
/*  Adapted from myprog.c
    Plots output text file via GNUPlot
*/
    char command[MAX_COMMAND_LEN];
    snprintf(command, sizeof(command), "%s %s", GNUPLOT_EXE, GNUPLOT_SCRIPT );
    system(command);
}

static double** Initialise2DArray(int length){
/*  Initialises 2d array with a row length of DIMS (n dimension vector) and a column length of 'length'
*/
    double** array = xmalloc(length*sizeof(double*));
    double* vals = xmalloc(DIMS*length*sizeof(double));
    for(int i=0;i<length;i++){      //start 2d array
        array[i] = vals + i*DIMS;
    }
    return array;
}

void Free2DArray(double** array, int bodyNo){
/*  Frees 2D array created by Initialise2DArray()
*/
    free(*array);
    free(array);
}

static double VectorMagnitude(double vec[DIMS]){
/*  Returns the magnitude of a DIMS long vector
*/
    double mag = 0;
    for(int i=0;i<DIMS;i++){
        mag += vec[i]*vec[i];
    }
    return sqrt(mag);
}

static double* VectorSubtract(double vec1[DIMS], double vec2[DIMS]){
/*  Returns the distance between two DIMS value arrays in the form of a DIMS long array
*/
    double* distanceVec = xmalloc(DIMS*sizeof(double));
    for(int i=0;i<DIMS;i++){
        distanceVec[i] = vec1[i]-vec2[i];
    }
    return distanceVec;
}

static double* VectorNormalise(double vec[DIMS]){
/*  Returns input vector scaled so its magnitude is one
*/
    double* vecN = xmalloc(DIMS*sizeof(double));
    double mag = VectorMagnitude(vec);
    for(int i=0;i<DIMS;i++){
        vecN[i] = vec[i]/mag;
    }
    return vecN;
}

static double KineticEnergy(BODY* bodies, int bodyNo){
/*  Returns total kinetic energy of BODY array
    Using equation kinetic energy E=1/2*m*v^2
*/
    double energy = 0;
    for(int i=0;i<bodyNo;i++){
        energy += 0.5*bodies[i].mass*pow(VectorMagnitude(bodies[i].vel),2);
    }
    return energy;
}

static double PotentialEnergy(BODY* bodies, int bodyNo){
/*  Returns total potential energy of BODY array
    Uses equation U = -G*M*m/r
*/
    double energy = 0;
    double* vecSub;
    for(int i=0;i<bodyNo;i++){
        for(int j=0;j<i;j++){ //So not repeating combinations of bodies
            if(i==j) fprintf(stderr,"ERROR\n");
            vecSub = VectorSubtract(bodies[i].pos,bodies[j].pos);
            energy += -BIG_G*bodies[i].mass*bodies[j].mass/VectorMagnitude(vecSub);
            free(vecSub);
        }
    }
    return energy;
}

static double TotalEnergy(BODY* bodies, int bodyNo){
/*  Returns total energy of the bodies system
*/
    return KineticEnergy(bodies,bodyNo)+PotentialEnergy(bodies,bodyNo);
}

static double** CalcPositions(BODY* bodies, int bodyNo, double timeStep){
/*  Returns 2D array created by 'Initialise2DArray' with values corresponding to the new positions of the bodies
    after 'timeStep' amount of time (Velocity Verlet method)
*/
    double** posF = Initialise2DArray(bodyNo);
    for(int i=0;i<bodyNo;i++){
        for(int j=0;j<DIMS;j++){
            posF[i][j] = bodies[i].pos[j] + bodies[i].vel[j]*timeStep + 0.5*bodies[i].acc[j]*pow(timeStep,2);
        }
    }
    return posF;
}

static double** CalcAccelerations(BODY* bodies, int bodyNo){
/*  Returns 2D array of new accelerations calculated via Velocity-Verlet method
    a = GM/r^2
*/
    double** accFinal = Initialise2DArray(bodyNo);
    double* r;
    double accMag = 0;
    double* vecNorm;
    for(int i=0;i<bodyNo;i++){
        for(int j=0;j<bodyNo;j++){
            if(i != j){
                r = VectorSubtract(bodies[i].pos, bodies[j].pos);
                accMag = -BIG_G * bodies[j].mass / (pow(VectorMagnitude(r),2)+SOFT_L*SOFT_L);
                vecNorm = VectorNormalise(r);
                for(int k=0;k<DIMS;k++){
                    accFinal[i][k] += accMag * vecNorm[k];
                }
                free(r);
                free(vecNorm);
            }
        }
    }
    return accFinal;
}

static double** CalcVelocities(BODY* bodies, int bodyNo, double** accNew, double timeStep){
/*  Returns 2D array of new velocities calculated via Velocity-Verlet method
*/
    double** velF = Initialise2DArray(bodyNo);
    for(int i=0;i<bodyNo;i++){
        for(int j=0;j<DIMS;j++){
            velF[i][j] = bodies[i].vel[j] + 0.5*timeStep*(bodies[i].acc[j]+accNew[i][j]);
        }
    }
    return velF;
}

static int Gfunc(double t, const double y[], double f[], void *params){
/*  Used in the GSL odeiv2 functions, Runge Kutta Method
    System of 2*DIMS differential equations
*/
    PARAMS p = *(PARAMS*)params;
    double constTerm;
    double yPos[DIMS];
    double posR[DIMS];

    for(int i=0;i<DIMS;i++){    //sets position of main body
        yPos[i] = y[i];
    }
    for(int i=DIMS;i<2*DIMS;i++){
        f[i] = 0;
    }

    for(int i=0;i<p.bodyNo;i++){    //Loops over each force applying body
        if(i != p.bodyID){
            double* vecSub = VectorSubtract(yPos,p.bodies[i].pos);
            for(int j=0;j<DIMS;j++){    // Sets relative position between objects
                posR[j] = vecSub[j];
            }
            constTerm = -(BIG_G*p.bodies[i].mass)/pow(posR[0]*posR[0]
                                                    + posR[1]*posR[1]
                                                    + posR[2]*posR[2] + SOFT_L*SOFT_L, 3.0/2.0);
            for(int j=0;j<DIMS;j++){
                f[j] = y[j+DIMS];
                f[j+DIMS] += posR[j]*constTerm;
            }
            free(vecSub);
        }
    }
    return 0;
}

static void VelocityVerletMethod(BODY* bodies, int bodyNo, double timeStep, int iterNo){
/*  Performs the Velocity Verlet method of solving the n body simulation
    Writes every FILE_PRINT_NO position of the bodies to file OUT_FILENAME
    Attempts to calculate the orbital periods for each body and prints these to stdout
*/
    FILE *file;
    file = fopen(OUT_FILENAME,"w");
    double** accNew = Initialise2DArray(bodyNo);    //1st time setup
    double** posNew = Initialise2DArray(bodyNo);
    double** velNew = Initialise2DArray(bodyNo);
    double pos_Zero[bodyNo];
    int halfway[bodyNo];    // If the period is past the half way point
    int recoveredP[bodyNo]; // If period for each body has been found

    for(int i=0;i<bodyNo;i++){  //Set up period recovery variables
        pos_Zero[i] = bodies[i].pos[1];
        halfway[i] = 0;
        recoveredP[i] = 0;
    }
    printf("Orbital Periods Retrieved (Days):\n");

    accNew = CalcAccelerations(bodies, bodyNo); //  Initial Calculation of accelerations
    for(int i=0;i<iterNo;i++){
        if(i){
            accNew = Initialise2DArray(bodyNo);
            posNew = Initialise2DArray(bodyNo);
            velNew = Initialise2DArray(bodyNo);
        }
        posNew = CalcPositions(bodies, bodyNo, timeStep);
        accNew = CalcAccelerations(bodies, bodyNo);
        velNew = CalcVelocities(bodies, bodyNo, accNew, timeStep);
        for(int j=0;j<bodyNo;j++){
            for(int k=0;k<DIMS;k++){
                bodies[j].pos[k] = posNew[j][k];    //assign new values to old values
                bodies[j].vel[k] = velNew[j][k];
                bodies[j].acc[k] = accNew[j][k];
                if(i%FILE_PRINT_NO == 0){
                    if(BODY_REF_FRAME == -1){
                        fprintf(file, "%lf\t", bodies[j].pos[k]);
                    }else{
                        fprintf(file, "%lf\t", bodies[j].pos[k] - bodies[BODY_REF_FRAME].pos[k]);
                    }
                }
                if(pos_Zero[j] >= bodies[j].pos[1] && i > 10 && halfway[j] == 0){   //Checks if body's y coord goes into -ve values
                    halfway[j] = 1;
                }
                if(pos_Zero[j] <= bodies[j].pos[1] && halfway[j] == 1 && recoveredP[j] == 0){   //Calculates period when y value becomes +ve again
                    printf("    %s:         %.2lf\n", bodies[j].name, timeStep*i/(60*60*24));
                    recoveredP[j] = 1;
                }
            }
            if(i%FILE_PRINT_NO == 0){fprintf(file,"\n");}
        }
        Free2DArray(posNew, bodyNo);
        Free2DArray(velNew, bodyNo);
        Free2DArray(accNew, bodyNo);
    }
    fclose(file);
}

static void RungeKutta4Method(BODY* bodies, int bodyNo, double timeStep, int iterNo){
/*  Performs the 4th order Runge Kutta method of solving the n body simulation
    Writes every FILE_PRINT_NO position of the bodies to file OUT_FILENAME
    Attempts to calculate the orbital periods for each body and prints these to stdout
*/
    double y[2*DIMS];
    double t = 0.0;
    int s = 0;
    int numSteps = 1;
    double pos_Zero[bodyNo];
    double period_time;     // Time of the period
    int halfway[bodyNo];    // If the period is past the half way point
    int recoveredP[bodyNo]; // If period for each body has been found
    PARAMS p;
    FILE *file;

    file = fopen(OUT_FILENAME,"w");

    for(int i=0;i<bodyNo;i++){
        pos_Zero[i] = bodies[i].pos[1];
        halfway[i] = 0;
        recoveredP[i] = 0;
    }
    printf("Orbital Periods Retrieved (Days):\n");

    p.bodies = bodies;
    p.bodyNo = bodyNo;

    for(int i=0;i<iterNo;i++){
        for(int j=0;j<bodyNo;j++){  // Current Body
            p.bodyID = j;
            gsl_odeiv2_system sys = {Gfunc, NULL, DIMS*2, &p};
            gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, 1e-3, 1e-8, 1e-8);
            for(int k=0;k<DIMS;k++){
                y[k] = bodies[j].pos[k];
                y[k+DIMS] = bodies[j].vel[k];
            }
            s = gsl_odeiv2_driver_apply_fixed_step(driver, &t, timeStep, numSteps, y);
            if(s != 0){
                fprintf(stderr,"ERROR: RK4 process returned %d\n", s);
                exit(EXIT_FAILURE);
            }
            for(int k=0;k<DIMS;k++){    // Update body pos/vel
                bodies[j].pos[k] = y[k];
                bodies[j].vel[k] = y[k+DIMS];
                if(i%FILE_PRINT_NO == 0){
                    if(BODY_REF_FRAME == -1){
                        fprintf(file, "%lf\t", bodies[j].pos[k]);
                    }else{
                        fprintf(file, "%lf\t", bodies[j].pos[k] - bodies[BODY_REF_FRAME].pos[k]);
                    }
                }
                if(pos_Zero[j] >= bodies[j].pos[1] && i > 10 && halfway[j] == 0){
                    halfway[j] = 1;
                }
                if(pos_Zero[j] <= bodies[j].pos[1] && halfway[j] == 1 && recoveredP[j] == 0){
                    period_time = timeStep*i/(60*60*24);
                    printf("    %s:         %.2lf\n", bodies[j].name, period_time);
                    recoveredP[j] = 1;
                }
            }
            if(i%FILE_PRINT_NO == 0){
                fprintf(file,"\n");
            }
            gsl_odeiv2_driver_free(driver);
        }
    }
    fclose(file);
}

int main(int argc, char **argv){
/*  Reads in command line arguments
    Performs error checking on input
    Prints calculated properties of the system
*/
    double timeStep = STD_STEP; //seconds
    double timeLen = STD_TIME_LEN; //seconds
    double startEnergy = 0;
    double endEnergy = 0;
    int iterNo = 1000000;
    int bodyNo = 0;
    int optionFlag = -1;
    int optionIndex = 0;
    BODY* bodies;
    char filename[MAX_LINE_SIZE];
    clock_t startT = clock();

    while(1){      //   This whole "while" loop is adapted from mat_gen.c
        static struct option long_options[] ={
            {"velocityVerlet", required_argument, 0, 'r'},
            {"rungeKutta", required_argument, 0, 'v'},
            {"timeStep", required_argument, 0, 't'},
            {"timeScale", required_argument, 0, 's'},
            {0, 0, 0, 0}
        };
        int c = getopt_long(argc, argv, "r:v:t:s:", long_options, &optionIndex);
        if (c == -1) break;
        switch (c){
            case 0:
                if (long_options[optionIndex].flag != 0)
                    break;
                fprintf(stderr,"option %s", long_options[optionIndex].name);
                if(optarg)
                    printf (" with arg %s", optarg);
                printf ("\n");
                break;
            case 'r':
                optionFlag = 1;
                printf("Method:     Runge Kutta 4th Order\n");
                snprintf(filename,sizeof(filename),"%s",optarg);
                break;
            case 'v':
                optionFlag = 0;
                printf("Method:     Velocity Verlet\n");
                snprintf(filename,sizeof(filename),"%s",optarg);
                break;
            case 't':
                timeStep = atof(optarg);
                if(timeStep <= 0){
                    fprintf(stderr, "Warning: Invalid time step value: %s, set to default of %.2lfs\n",optarg,STD_STEP);
                    timeStep = STD_STEP;
                }
                break;
            case 's':
                timeLen = atof(optarg);
                if(timeLen <= 0){
                    fprintf(stderr, "Warning: Invalid time scale value: %s, set to default of %.1Es\n",optarg,STD_TIME_LEN);
                    timeLen = STD_TIME_LEN;
                }
                break;
            case ':':
                printf("ERROR: No argument specified\n");
                break;
            default:
                printf("ERROR: Invalid option\n");
        }
    }
    if (optind < argc) {    // Checks for an invalid argument
        fprintf(stderr,"Invalid Argument: ");
        while (optind < argc) {
            fprintf(stderr, "%s ", argv[optind++]);
        }
        fprintf(stderr, "\n");
    }

    iterNo = (int)floor(timeLen/timeStep);
    if(iterNo <= 0){    // Checks for valid iteration number
        fprintf(stderr,"ERROR: Invalid selection of time step and time scale");
        exit(EXIT_FAILURE);
    }

    bodyNo = ReadFileBodyNo(filename);  // Reads in file data
    if(bodyNo == -1){       //  Checks for valid number of bodies
        fprintf(stderr,"ERROR: Invalid File\n");
        exit(EXIT_FAILURE);
    }
    else if(bodyNo < 2){
        fprintf(stderr,"ERROR: At least two valid bodies needed in file\n");
        exit(EXIT_FAILURE);
    }
    bodies = ReadFileData(filename, bodyNo);
    printf("Number of bodies detected:  %i\n    (",bodyNo);
    for(int i=0;i<bodyNo-1;i++){
        printf("%s, ",bodies[i].name);
    }
    printf("%s)\n", bodies[bodyNo-1].name);
    startEnergy = TotalEnergy(bodies, bodyNo);

    if(optionFlag == -1){   // If no methods are selected
        fprintf(stderr,"ERROR: No method option selected\n");
        exit(EXIT_FAILURE);
    }
    else if(optionFlag == 0){
        VelocityVerletMethod(bodies, bodyNo, timeStep, iterNo);
        GnuPlot();
    }
    else if(optionFlag == 1){
        RungeKutta4Method(bodies, bodyNo, timeStep, iterNo);
        GnuPlot();
    }
    //Print out Energy properties
    endEnergy = TotalEnergy(bodies, bodyNo);
    printf("Total Energy (J):\n");
    printf("    Initial:         %.3E\n",startEnergy);
    printf("    Final:           %.3E\n",endEnergy);
    printf("    Change:          %.3E\n",endEnergy-startEnergy);
    printf("    Percent Error:   %.3lg%%\n",100*fabs((endEnergy-startEnergy)/startEnergy));

    free(bodies);
    clock_t endT = clock();
    printf("Execution Time (s):  %.4lfs\n",(double)(endT-startT)/CLOCKS_PER_SEC);
    return 0;
}



