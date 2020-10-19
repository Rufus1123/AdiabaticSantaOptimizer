// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;
using Microsoft.Quantum.Simulation.Simulators;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace Microsoft.Quantum.Samples.Ising
{
    class Program
    {
        static void Main(string[] args)
        {
            #region Basic Definitions

            // We start by loading the simulator that we will use to run our Q# operations.
            var qsim = new QuantumSimulator();

            // In all the following, we use this coefficient for coupling to the transverse
            // field.
            var hxCoeff = 1.0;

            // As we are using a Trotter–Suzuki decomposition as our simulation algorithm,
            // we will need to pick a timestep for the simulation, and the order of the
            // integrator. The optimal timestep needs to be determined empirically, and
            // we find that the following choice works well enough.
            var trotterStepSize = 0.05;
            var trotterOrder = 2;

            // Let us now simulate time-evolution by interpolating between the initial
            // Hamiltonian with the |+〉 product state as the ground state, and the target
            // Hamiltonian. For the uniform Ising model, the ground state of the target
            // Hamiltonian should have all spins pointing in the same direction. If we
            // interpolate between these Hamiltonians slowly enough, the initial ground
            // state will continuously deform into the ground state of the target
            // Hamiltonian

            // Let now interpolate between these Hamiltonians more slowly.
            var adiabaticTime = 500.0;

            #endregion

            #region Ising model simulations

            var penalty = 20.0;
            var segmentCosts = new double[] { 4.70, 9.09, 9.03, 5.70, 8.02, 1.71 };
            var numQubits = 6;

            // The phase estimation algorithm requires us to choose the duration of time-
            // evolution in the oracle it calls, and the bits of precision to which we
            // estimate the phase. Note that the error of the energy estimate is typically
            // rescaled by 1 / `qpeStepSize`
            var qpeStepSize = 0.1;
            var nBitsPrecision = 10;

            var linearTerm = segmentCosts.Select(c => 4.0 * penalty - 0.5 * c).ToArray();
            var quadraticTerm = new double[numQubits * numQubits].Select((c, i) => i == 2 || i == 9 || i == 29 ? penalty : 2.0 * penalty).ToArray();

            var hZFinal = new QArray<Double>(linearTerm);
            var jZZFinal = new QArray<Double>(quadraticTerm);
            CalculateCache(quadraticTerm, linearTerm, segmentCosts);

            Console.WriteLine("\nIsing model parameters:");
            Console.WriteLine(
                $"\t{numQubits} sites\n" +
                $"\t{hxCoeff} transverse field coefficient\n" +
                $"\t{adiabaticTime} time-interval of interpolation\n" +
                $"\t{trotterStepSize} simulation time step \n" +
                $"\t{trotterOrder} order of integrator\n");

            // Let us use this opportunity to test the adiabatic evolution as written using
            // more library functions.
            for (int rep = 0; rep < 10; rep++)
            {
                var data = GeneralIsingAdiabaticAndMeasureBuiltIn.Run(qsim, 
                    numQubits, hxCoeff, hZFinal, jZZFinal, 
                    adiabaticTime, trotterStepSize, trotterOrder, 
                    qpeStepSize, nBitsPrecision).Result;

                var measuredState = data.ToArray().Select(x => x.ToString() == "Zero" ? "0" : "1");

                Console.Write($"State: {string.Join("", measuredState)} \n");

                var expectationValues = CalculateExpectationValuesFromFile("..\\..\\..\\QubitsDump.txt", numQubits);

                Console.WriteLine($"Energy after evolution: {expectationValues.Item1}\nExpected Costs: {expectationValues.Item2}\nCosts for measured state: {CalculateCostForStateCache(measuredState.Select(q => q.Equals("1")).ToArray(), segmentCosts)}");
            }

            Console.WriteLine("Press Enter to continue...");
            Console.ReadLine();

            #endregion
        }




        // DEBUG
        public static Dictionary<string, double> energyStateCache = new Dictionary<string, double>();
        public static Dictionary<string, double> costStateCache = new Dictionary<string, double>();

        public static void CalculateCache(double[] J, double[] h, double[] segmentCosts)
        {
            var numQubits = h.Length;
            for (int i = 0; i < 1 << numQubits; i++)
            {
                var key = string.Join("", System.Convert.ToString(i, 2).PadLeft(numQubits, '0'));
                var qubits = key.Select(s => s.Equals('1')).ToArray();

                var energyForState = CalculateEnergyForStateCache(qubits, J, h);
                energyStateCache.Add(key, energyForState);

                var costForState = CalculateCostForStateCache(qubits, segmentCosts);
                costStateCache.Add(key, costForState);
            }
        }

        private static double CalculateEnergyForStateCache(bool[] qubits, double[] J, double[] h)
        {
            var numQubits = qubits.Length;
            var z = new double[numQubits];

            for (var i = 0; i < numQubits; i++)
            {
                z[i] = qubits[i] ? 1 : -1;
            }

            var totalEnergy = 0.0;

            for (var i = 0; i < numQubits; i++)
            {
                for (var j = i + 1; j < numQubits; j++)
                {
                    totalEnergy += J[i*numQubits + j] * z[i] * z[j];
                }

                totalEnergy += h[i] * z[i];
            }

            return totalEnergy;
        }

        private static double CalculateCostForStateCache(bool[] qubits, double[] segmentCosts)
        {
            var numQubits = qubits.Length;
            var z = new double[numQubits];

            for (var i = 0; i < numQubits; i++)
            {
                z[i] = qubits[i] ? 0 : 1;
            }

            var totalCost = 0.0;

            for (var i = 0; i < numQubits; i++)
            {

                totalCost += segmentCosts[i] * z[i];
            }

            return totalCost;
        }

        public static (double, double) CalculateExpectationValuesFromFile(string fileName, int numQubits)
        {
            string[] lines = File.ReadAllLines(fileName);

            var expectedEnergy = 0.0;
            var expectedCost = 0.0;

            foreach (var (line, i) in lines.Select((value, index) => (value, index - 1)))
            {
                var lineArr = line.Split('\t');
                if (lineArr.Length > 1)
                {

                    var realPart = double.Parse(new Regex(@"-?[0-9][0-9,\.]+").Match(lineArr[1]).Value);
                    var imaginaryPart = double.Parse(new Regex(@"-?[0-9][0-9,\.]+\si").Match(lineArr[1]).Value.Replace(" i", ""));
                    var probability = realPart * realPart + imaginaryPart * imaginaryPart;

                    string key = string.Join("", System.Convert.ToString(i, 2).PadLeft(numQubits, '0'));
                    var energy = energyStateCache[key];
                    expectedEnergy += probability * energy;

                    var cost = costStateCache[key];
                    expectedCost += probability * cost;
                }
            }

            return (expectedEnergy, expectedCost);
        }
    }
}
