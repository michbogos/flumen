#ifndef FLUMEN
#define FLUMEN

#define GRIDX 128
#define GRIDY 128
#define G -9.81
#define RHO 1

//Horizontal at center of vertical
float u[GRIDX][GRIDY];
//Vertical at center of horizontal
float v[GRIDX][GRIDY];

//Pressure
float p[GRIDX][GRIDY];

void update(float** u, float** v, float dt){
    for(int i = 0; i < GRIDY;i++){
        for(int j = 0; j < GRIDX; j++){
            //Velocity
            v[i][j] += G*dt;
        }
    }
}

void project(float** u, float** v, float** p, float dt){
    for(int n = 0; n < 20; n++){
        for(int i = 0; i < GRIDY; i++){
            for(int j = 0; j < GRIDX; j++){
                if(i > 0 && i < GRIDY-1){
                    if(j > 0 && j < GRIDX-1){
                        float d = u[i+1][j]-u[i][j]+v[i][j+1]-v[i][j];
                        u[i][j]+=d/4;
                        u[i+1][j] -= d/4;
                        v[i][j]+=d/4;
                        v[i][j+1] -= d/4;
                        p[i][j] += d/4*RHO*1/GRIDX/dt;
                    }
                }
            }
        }
    }
}

void advect(float** u, float** v, float dt){
    float vbar = 
}



#endif