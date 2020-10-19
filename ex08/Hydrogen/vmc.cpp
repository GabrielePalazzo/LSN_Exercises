#include <iostream>
#include "vmc.h"

using namespace std;

void Minimize(WaveFunction* wf, double startPosition, Random& generator, double* gradient, int step, double* parametersBest, double& energyBest) {
    
    //double temperature = ((500 - step/2) * (500 - step/2) * (500 - step/2))/1000000000.0 + 1/2000.0;
    double temperature = exp(-step/(double)80);
    
    int accepted = 0;
    for (int i = 0; i < equilibriumSteps; i++) {
        Metropolis_Step(startPosition, generator, samplingType, wf, accepted);
    }
    Simulate_Annealing(wf, startPosition, generator, gradient, step, parametersBest,  energyBest, temperature);
    /*
    int accepted = 0;
    for (int i = 0; i < equilibriumSteps; i++) {
        Metropolis_Step(startPosition, generator, samplingType, wf, accepted);
    }
    Gradient_Descent(wf, startPosition, generator, gradient, step, parametersBest,  energyBest);
    */
}

void Simulate_Annealing(WaveFunction* wf, double startPosition, Random& generator, double* gradient, int step, double* parametersBest, double& energyBest, double temperature) {
    double oldEnergy = 0, newEnergy = 0;
    int accepted = 0;
    int stepsPerBlock = numberOfSteps / numberOfBlocks;
    double averageEnergy = 0;
    double averageEnergySquared = 0;
    double x = startPosition;
    for (int i = 0; i < numberOfBlocks; i++) {
        double energyInsideBlock = 0;
        for (int j = 0; j < stepsPerBlock; j++) {
            Metropolis_Step(x, generator, samplingType, wf, accepted);
            energyInsideBlock += wf->Energy(x);
        }
        energyInsideBlock /= stepsPerBlock;
        averageEnergy += energyInsideBlock;
        averageEnergySquared += energyInsideBlock * energyInsideBlock;
        
    }
    oldEnergy = averageEnergy/numberOfBlocks;
    x = startPosition;
    averageEnergy = 0;
    averageEnergySquared = 0;
    
    
    double oldMu = wf->Get_mu();
    double oldSigma = wf->Get_sigma();
    
    if (oldEnergy < energyBest) {
        parametersBest[0] = oldMu;
        parametersBest[1] = oldSigma;
        energyBest = oldEnergy;
    }
    
    double newMu = 0;
    double newSigma = 0;
    
    //if (step < 1) {
    gradient[0] = (generator.Rannyu() - 0.5) * 1;
    gradient[1] = (generator.Rannyu() - 0.5) * 0.1;
    //}
    /*
     else {
     gradient[0] /= step;
     }*/
    //cout << endl << gradient[0] << endl;
    if (step < 1) {
        newMu = oldMu + gradient[0];
        newSigma = oldSigma + gradient[1];
    }
    else {
        newMu = oldMu + gradient[0];
        newSigma = oldSigma + gradient[1];
    }
    
    wf->Set_mu(newMu);
    wf->Set_sigma(newSigma);
    
    for (int i = 0; i < numberOfBlocks; i++) {
        double energyInsideBlock = 0;
        for (int j = 0; j < stepsPerBlock; j++) {
            Metropolis_Step(x, generator, samplingType, wf, accepted);
            energyInsideBlock += wf->Energy(x);
        }
        energyInsideBlock /= stepsPerBlock;
        averageEnergy += energyInsideBlock;
        averageEnergySquared += energyInsideBlock * energyInsideBlock;
        
    }
    newEnergy = averageEnergy/numberOfBlocks;
    
    if (newEnergy < energyBest) {
        parametersBest[0] = newMu;
        parametersBest[1] = newSigma;
        energyBest = newEnergy;
    }
    
    wf->Set_mu(oldMu);
    wf->Set_sigma(oldSigma);
    
    if (generator.Rannyu() < exp((oldEnergy - newEnergy)/temperature)) {
        wf->Set_mu(newMu);
        wf->Set_sigma(newSigma);
        
        //cout << temperature << " " << newMu << " " << newSigma << " " << newEnergy << endl;
    }
}

void Gradient_Descent(WaveFunction* wf, double startPosition, Random& generator, double* gradient, int step, double* parametersBest, double& energyBest) {
    double oldEnergy = 0;
    int accepted = 0;
    int stepsPerBlock = numberOfSteps / numberOfBlocks;
    double averageEnergy = 0;
    double averageEnergySquared = 0;
    double averageMu = 0;
    double averageMuSquared = 0;
    double averageSigma = 0;
    double averageSigmaSquared = 0;
    double x = startPosition;
    double gradMu = 0, gradSigma = 0;
    for (int i = 0; i < numberOfBlocks; i++) {
        double energyInsideBlock = 0;
        double muInsideBlock = 0;
        double sigmaInsideBlock = 0;
        for (int j = 0; j < stepsPerBlock; j++) {
            Metropolis_Step(x, generator, samplingType, wf, accepted);
            energyInsideBlock += wf->Energy(x);
            muInsideBlock += wf->Grad_mu(x);
            sigmaInsideBlock += wf->Grad_sigma(x);
        }
        energyInsideBlock /= stepsPerBlock;
        muInsideBlock /= stepsPerBlock;
        sigmaInsideBlock /= stepsPerBlock;
        
        averageEnergy += energyInsideBlock;
        averageEnergySquared += energyInsideBlock * energyInsideBlock;
        averageMu += muInsideBlock;
        averageMuSquared += muInsideBlock * muInsideBlock;
        averageSigma += sigmaInsideBlock;
        averageSigmaSquared += sigmaInsideBlock * sigmaInsideBlock;
        
    }
    oldEnergy = averageEnergy/numberOfBlocks;
    gradMu = averageMu/numberOfBlocks;
    gradSigma = averageSigma/numberOfBlocks;
    
    double oldMu = wf->Get_mu();
    double oldSigma = wf->Get_sigma();
    
    if (oldEnergy < energyBest) {
        parametersBest[0] = oldMu;
        parametersBest[1] = oldSigma;
        energyBest = oldEnergy;
    }
    
    double alpha = 0.3;
    double beta = 0.9;
    
    double newMu = 0;
    double newSigma = 0;/*
    if (step < 1) {
        gradient[0] = (generator.Rannyu() - 0.5) * 2.0;
        //gradient[0] = beta * gradMu;
        gradient[0] = 1;
        gradient[1] = (generator.Rannyu() - 0.5) * 2.0;
        //gradient[1] = beta * gradSigma;
        gradient[1] = 1;
    }
    else {
        //gradient[0] /= step;
        //gradient[1] /= step;
    }*/
    
    
    double steepness = gradMu;
    int sign = 1;
    if (steepness < 0) sign = -1;
    gradient[0] = beta * gradient[0] + sign / (1 + sign * steepness);
    //gradient[0] = beta * gradient[0] + steepness;
    
    newMu = oldMu - alpha * gradient[0];
    //newMu = newMu - sign * newMu / 2.0;
    
    
    steepness = gradSigma;
    sign = 1;
    if (steepness < 0) sign = -1;
    
    gradient[1] = beta * gradient[1] + sign / (1 + sign * steepness);
    //gradient[1] = beta * gradient[1] + steepness;
    
    newSigma = oldSigma - alpha * gradient[1];
    double threshold = 0.001;
    if (gradMu < threshold && gradSigma < threshold) cout << "Ciao" << endl;
    else {
    wf->Set_mu(newMu);
    wf->Set_mu(newSigma);
    }
    //newSigma = oldSigma;
    /*
    wf->Set_mu(newMu);
    averageEnergy = 0;
    averageEnergySquared = 0;
    x = startPosition;
    for (int i = 0; i < numberOfBlocks; i++) {
        double energyInsideBlock = 0;
        for (int j = 0; j < stepsPerBlock; j++) {
            Metropolis_Step(x, generator, samplingType, wf, accepted);
            energyInsideBlock += wf->Energy(x);
        }
        energyInsideBlock /= stepsPerBlock;
        averageEnergy += energyInsideBlock;
        averageEnergySquared += energyInsideBlock * energyInsideBlock;
        
    }
    newEnergy = averageEnergy/numberOfBlocks;
    
    double steepness = ((newEnergy - oldEnergy) / (newMu - oldMu));
    int sign = 1;
    if (steepness < 0) sign = -1;
    gradient[0] = 0.5 * gradient[0] - 0.4 * sign * 1/(1 + sign * steepness) + (generator.Rannyu() - 0.5) * 0.2;
    //gradient[0] = 0.0 * gradient[0] - 0.01 * ((newEnergy - oldEnergy) / (newMu - oldMu));
    wf->Set_mu(oldMu);
    
    
    wf->Set_mu(newSigma);
    averageEnergy = 0;
    averageEnergySquared = 0;
    x = startPosition;
    for (int i = 0; i < numberOfBlocks; i++) {
        double energyInsideBlock = 0;
        for (int j = 0; j < stepsPerBlock; j++) {
            Metropolis_Step(x, generator, samplingType, wf, accepted);
            energyInsideBlock += wf->Energy(x);
        }
        energyInsideBlock /= stepsPerBlock;
        averageEnergy += energyInsideBlock;
        averageEnergySquared += energyInsideBlock * energyInsideBlock;
        
    }
    newEnergy = averageEnergy/numberOfBlocks;
    
    
    steepness = ((newEnergy - oldEnergy) / (newSigma - oldSigma));
    sign = 1;
    if (steepness < 0) sign = -1;
    
    gradient[1] = 0.5 * gradient[1] - 0.4 * sign * 1/(1 + sign * steepness) + (generator.Rannyu() - 0.5) * 0.2;
    //gradient[1] = 0.0 * gradient[1] + 0.01 * ((newEnergy - oldEnergy) / (newSigma - oldSigma));
    
    //wf->Set_mu(oldSigma);
    wf->Set_mu(newMu);
     */
    cout << "\t" << oldEnergy << endl;
    cout << newMu << " " << newSigma << endl;
    cout << "\t\t" << gradMu << " " << gradSigma << endl;
    cout << "\t\t" << gradient[0] << " " << gradient[1] << endl;
}


