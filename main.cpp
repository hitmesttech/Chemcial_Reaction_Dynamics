#include <cstdio>
#include <cmath>
#include <cstdlib>

using namespace std;

#define MAX_REACTION_CONSTANT 1024

int sample_frequency;
long double simu_step,simu_time;

int reactant_number;
int reaction_number;
int reactant_serial[MAX_REACTION_CONSTANT][MAX_REACTION_CONSTANT];
long double reaction_constant[MAX_REACTION_CONSTANT];
long double reaction_active_energy[MAX_REACTION_CONSTANT];
int reactant_involved[MAX_REACTION_CONSTANT];//每个反应都包含了什么多少种反应物
long double temperature;
long double chemometric[MAX_REACTION_CONSTANT][MAX_REACTION_CONSTANT];
long double dynamic_ratio[MAX_REACTION_CONSTANT][MAX_REACTION_CONSTANT];

char reactant_name[MAX_REACTION_CONSTANT][MAX_REACTION_CONSTANT];
//this function will cause bug. Use dynamic_list_alloc instead;

inline int getsite(char c,char *s,int j, FILE *in)
{
     while(fgets(s,j,in))
        {
            if(s[0]==c)
                break;
        }
}

inline void** dynamic_matrix_alloc(int a,int b,size_t j)
{
    int i;
    void **pt,*x;
    pt=new void*[a];
    for(i=0;i<a;i++)
    {
        (pt[i])=malloc(b*j);
    }
    return pt;
}


inline void** dynamic_list_alloc(int a,int b[],size_t j)
{
    int i;
    void **pt;
    pt=new void*[a];
    for(i=0;i<a;i++)
    {
        (pt[i])=malloc(b[i]*j);
    }
    return pt;
}
//this function calculate the sigma of one reaction
inline long double sigma_one_reaction_one_step(int reaction_serial,long double reactant_amount[])
{
    int i,j=reaction_serial;//j stands for the reaction serial;
    long double sigma=reaction_constant[reaction_serial]*
    exp(-reaction_active_energy[reaction_serial]*1000/(8.314*temperature) );
    for(i=0;i<reactant_involved[j];i++)
    {
        sigma *= pow(reactant_amount[ (reactant_serial[reaction_serial][i]) ],
        dynamic_ratio[reaction_serial][ (reactant_serial[reaction_serial][i]) ]);
    }
    return sigma;
}

inline int all_reaction_one_step(long double reactant_amount[],long double reactant_amount_delta[])
{
    int i,j;
    long double sigma;
    for(i=0;i<reactant_number;i++)
    {
        reactant_amount_delta[i]=0;
    }
    for(i=0;i<reaction_number;i++)
    {
        sigma=sigma_one_reaction_one_step(i,reactant_amount);
        for(j=0;j<reactant_involved[i];j++)
        {
            reactant_amount_delta[ (reactant_serial[i][j]) ]+=sigma*chemometric[i][ (reactant_serial[i][j]) ];
        }
    }
    for(i=0;i<reactant_number;i++)
    {
        reactant_amount[i]+=reactant_amount_delta[i];
    }
    
    return 0;
}

int simulation(FILE *ou,long double reactant_amount[])
{
    int j;
    long long i,tstep,tstep_higher;
    tstep=simu_time/simu_step;
    long double reactant_amount_delta[MAX_REACTION_CONSTANT];
    
    fprintf(ou,"Step\tTime");
    for(j=0;j<reactant_number;j++)
    {
        fprintf(ou,"\t%s",reactant_name[j]);
    }
    fprintf(ou,"\n");
    
    for(i=0;i<tstep;i++)
    {
        all_reaction_one_step(reactant_amount,reactant_amount_delta);
        if(!(i%sample_frequency))
        {
            printf("%% %lf completed, %lld step finished\r",( (double)i )/( (double)tstep ),tstep);
            
            fprintf(ou,"%lld\t%lf",i,i*simu_step);
            for(j=0;j<reactant_number;j++)
            {
                fprintf(ou,"\t%llf",reactant_amount[i]);
            }
            fprintf(ou,"'\n");
        }
    }
    
    return 0;
}

int simulation_today(FILE *ou[5],long double reactant_amount_inital[])
{
    int j,k;
    long long i,tstep,tstep_higher;
    tstep=simu_time/simu_step;
    long double reactant_amount[5][MAX_REACTION_CONSTANT],reactant_amount_delta[5][MAX_REACTION_CONSTANT];
    long double delta;
    
    for(j=0;j<reactant_number;j++)
    {
        reactant_amount[0][j]=reactant_amount[1][j]=reactant_amount[2][j]=reactant_amount[3][j]=reactant_amount[4][j]=reactant_amount_inital[j];
    }
    
    fprintf(ou[i],"Step\tTime");
    for(j=0;j<reactant_number;j++)
    {
        fprintf(ou[i],"\t%s",reactant_name[j]);
    }
    fprintf(ou[i],"\n");
    
    for(i=0;i<tstep;i++)
    {
            for(k=0;k<5;k++)
            {
                all_reaction_one_step(reactant_amount[k],reactant_amount_delta[0]);
                if(!(i%sample_frequency))
                {
                    printf("thread %d %% %lf completed, %lld step finished\r",k,( (double)i )/( (double)tstep ),tstep);

                    fprintf(ou[k],"%lld\t%lf",i,i*simu_step);
                    for(j=0;j<reactant_number;j++)
                    {
                        fprintf(ou[k],"\t%llf",reactant_amount[i][k]);
                    }
                    fprintf(ou[k],"\n");
                }
                for(j=0,delta=0;j<reactant_number;j++)
                {
                    delta+=(                     
                    fabs(reactant_amount[1][j]-reactant_amount[0][j])+
                    fabs(reactant_amount[2][j]-reactant_amount[0][j])+
                    fabs(reactant_amount[3][j]-reactant_amount[0][j])+
                    fabs(reactant_amount[4][j]-reactant_amount[0][j]) 
                    )/(5*reactant_amount[0][j]);
                }
                if(delta>1e-4)
                    printf("error in step %lld\n",i);
            }
    }
    
    return 0;
}

int getreactants(char *reactants_path)
{
    int i;
    char str_data[2048];
    FILE *in=fopen(reactants_path,"r");
    while(fgets(str_data,2048,in))
    {
        if(str_data[0]=='$')
            sscanf(str_data,"$reactant_number=",&reactant_number);
        else if(str_data[0]='//')
            continue;
        else 
        {
            sscanf(str_data,"%ld",&i);
            sscanf(str_data,"%ld%s",&i,reactant_name[i]);
        }
    } 
}

int getreactions(char reactions_path[])
{
    int i,reaction_serial,j;
    char str_data[2048];
    FILE *in=fopen(reactions_path,"r");
    fscanf(in,"$reaction_number=%ld",&reaction_number);
    for(i=0;i<reaction_number;i++)
    {
        
        getsite('&',str_data,2048,in);
        
        sscanf(str_data,"&&reaction%ld",&reaction_serial);
        
        getsite('$',str_data,2048,in);
        
        sscanf(str_data,"$activation_energy=%lf",reaction_active_energy+reaction_serial);
        
        getsite('$',str_data,2048,in);
        
        sscanf(str_data,"$reaction_constant=%lf",reaction_constant+reaction_serial);
        
        getsite('$',str_data,2048,in);
        
        sscanf(str_data,"$reactants_involved=%d",reactant_involved+reaction_serial);
        
        getsite('$',str_data,2048,in);
        
        for(j=0;j<reactant_involved[i];j++)
        {
            fscanf(in,"%ld",&reactant_serial[reaction_serial][j]);
        }
        
        getsite('$',str_data,2048,in);
        
        for(j=0;j<reactant_involved[i];j++)
        {
            fscanf(in,"%ld",&chemometric[reaction_serial][j]);
        }
        
        getsite('$',str_data,2048,in);
        
        for(j=0;j<reactant_involved[i];j++)
        {
            fscanf(in,"%ld",dynamic_ratio[reaction_serial][j]);
        }
    } 
}


int getcondition(char *reactants_path,long double reactant_amount[])
{
    int i;
    char str_data[2048];
    FILE *in=fopen(reactants_path,"r");
    while(fgets(str_data,2048,in))
    {
        if(str_data[0]=='$'&&str_data[1]=='r')
            sscanf(str_data,"reactant_number=",&reactant_number);
        else if(str_data[0]='//')
            continue;
        else 
        {
            sscanf(str_data,"%ld",&i);
            sscanf(str_data,"%ld%s",&i,reactant_amount+i);
        }
    } 
}

int main(int argc, char **argv)
{
    char reactants_path[2048];
    char reaction_path[2048];
    char initial_condition_path[2048];
    long double reactant_amount[MAX_REACTION_CONSTANT];
    int i;
    printf("Molecular Dynamics Simulation Set:\nChemical Reaction Dynamics Simulation Tools v2.0\n");
    printf("Developed by Teren Liu@Harbin Institute of Technology.\nPublished Under BSDv4 License\n");
    printf("Last updated:2015-12-9@Dongbo Wang's team, Department of Information Material Science,HIT");
    printf("First version published :2014-9-30.thanks to the HIT-iGEM team\n");
    
    printf("\nPlease input the reactant path:");
    gets(reactants_path);
    
    printf("\n Please input the reaction path:");
    gets(reaction_path);
    
    return 0;
}
