/*
Author : @Rajat Shah
Task: Solving N Queen using Genetic Algorithm.
Comments: Hard coded for 8x8 chess, can be modified for arbitrary size.
*/

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <sstream>
#include <cstring>
#include <map>
#include <list>
#include <queue>
#include <ctime>

using namespace std;

typedef struct
{
    string arrangement;
    int cost;
} individual;

typedef vector<individual*> population_type;

population_type population;
int chessBoardSize;
int initialPopulationCount = 10;

int fitnessValue(string arrangement)
{
    //initialize to a solution
    int fitness = (chessBoardSize*(chessBoardSize-1)) / 2;

    //removing pairs that lie on the same row and on the same diagonal
    for(int i=0; i<chessBoardSize; i++)
        for(int j=i+1; j<chessBoardSize; j++)
            if((arrangement[i] == arrangement[j]) || 
                (i-arrangement[i] == j-arrangement[j]) || 
                (i+arrangement[i] == j+arrangement[j])) 
            {
                fitness--;
            }

    return fitness;
}

individual* createNode()
{
    individual *newNode  = new individual;
    return newNode;
}

void generatePopulation()
{
    string sampleArrangement="";
    individual *temp;

    for(int i=1; i<=chessBoardSize; i++)
    {
        ostringstream ostr;
        ostr<<i;
        sampleArrangement += ostr.str();
    }

    //adds entries to population list
    for(int i=0; i<initialPopulationCount; i++)
    {
        random_shuffle( sampleArrangement.begin(), sampleArrangement.end());
        temp   = createNode();
        temp->arrangement = sampleArrangement;
        temp->cost = fitnessValue(sampleArrangement);
        population.push_back(temp);
    }
}

individual* reproduce(individual *x, individual *y)
{
    individual *child = createNode();
    int n = chessBoardSize;
    int c = rand() % n;
    child->arrangement = (x->arrangement).substr(0,c) + (y->arrangement).substr(c,n-c+1);
    child->cost = fitnessValue(child->arrangement);
    return child;
}

individual* mutate(individual *child)
{
    int randomQueen = rand()%(chessBoardSize)+1;
    int randomPosition= rand()%(chessBoardSize)+1;
    child->arrangement[randomQueen] = randomPosition+48;
    return child;
}

int randomSelection()
{
    int randomPos = rand()%population.size() %2;
    return randomPos;
}

bool isFit(individual *test)
{
    if(fitnessValue(test->arrangement)==((chessBoardSize*(chessBoardSize-1))/2))
        return true;
    return false;
}

bool comp(individual *a, individual*b)
{
    return(a->cost > b->cost);
}

individual* GA()
{
    int randomNum1,randomNum2;
    individual *individualX,*individualY,*child;
    bool found =0;
    while(!found)
    {
        population_type new_population;
        for(unsigned int i=0; i<population.size(); i++)
        {
            sort(population.begin(),population.end(),comp);

            randomNum1 = randomSelection();
            individualX = population[randomNum1];

            randomNum2 =randomSelection();
            individualY = population[randomNum2];

            child = reproduce(individualX,individualY);

            if(rand() % 2 == 0)     //random probability
                child = mutate(child);

            if(isFit(child))
            {
                found=1;
                return child;
            }
            new_population.push_back(child);
        }
        population = new_population;
    }
    return child;
}

void initialize()
{
    srand(time(0));     //to ensure perfect randomness
    chessBoardSize=15;
}

int main()
{
    initialize();
    clock_t start_time, end_time;           //to keep a track of the time spent in computing

    map<string,int> solutionsFound;
    int maxSolutions = 92, numFound = 0;

    start_time = clock();
    cout << "*Returns the column number corresponding to the row at the index*" << endl << endl;

    while(numFound!=maxSolutions)
    {
        generatePopulation();
        individual *solution = GA();
        if(!solutionsFound[solution->arrangement])
        {
            solutionsFound[solution->arrangement] = 1;
            cout<<"Possible Solution #"<<(++numFound)<<":\t"<<solution->arrangement<<endl;
        }
    }

    end_time = clock();
    cout << "Time required for execution: \t"  << 1000*((double)(end_time-start_time)/CLOCKS_PER_SEC) << " milliseconds." << "\n\n";
    return 0;
}
