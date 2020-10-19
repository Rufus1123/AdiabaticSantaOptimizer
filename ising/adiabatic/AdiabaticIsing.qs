// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Samples.Ising {
    open Microsoft.Quantum.Characterization;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Simulation;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Measurement;

    //////////////////////////////////////////////////////////////////////////
    // Introduction //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    // In this sample, the generator representation of the Ising model
    // that we constructed in the Ising Generators Sample will be used as the
    // input to simulation algorithms in the canon. We will use these
    // simulation algorithms to realize adiabatic state preparation.

    // In adiabatic state preparation, we interpolate between two Hamiltonians.
    // We begin with the initial Hamiltonian Hᵢ, which has an easy-to-prepare
    // ground state |ψᵢ〉. This Hamiltonian is then continuously deformed into
    // the target Hamiltonian Hₜ, with the desired ground state |ψₜ〉 that is
    // typically more difficult to prepare.

    // These define the interpolated Hamiltonian

    // H(s) = (1-s) Hᵢ + s Hₜ,

    // where s ∈ [0,1] is a schedule parameter. Typically, the schedule
    // parameter is linearly related to the physical time t ∈ [0,T]. For
    // instance, if interpolation between the Hamiltonians occur over a
    // total time T, one may define

    // s = t / T

    // This is of course not the only possible choice, and may be generalized
    // to s = f(t), where f is some arbitrary function that satisfies f(0) = 0
    // and f(T) = 1. Crucially, one necessary condition of the procedure is
    // that H(s) is continuous with respect to s, and thus f(t) is a
    // continuous function.

    // By performing time-evolution by H(s) while slowly varying the schedule
    // from 0 to 1 over physical time T, the initial ground state |ψᵢ〉 remains
    // an instantaneous ground state of H(s), and when s = 1, is then
    // transformed into the target ground state |ψₜ〉. The probability of
    // success improves the larger T is, equivalently, the more
    // slowly s is varied per unit of physical time.

    // In many cases, the optimal physical time of the interpolation T is
    // determined empirically. Though the worst-case rate of varying s can be
    // obtained from the gap of the Hamiltonian, which is the difference in
    // energy between the instantaneous ground state and the first excited
    // state, computing the gap is in general an intractable problem.

    // In other situations, one may also choose a non-linear schedule f(t)
    // which may impart desirable properties, such as reduced error in state
    // preparation, or even allow for shorter T. Choosing the optimal f is,
    // however, a very difficult problem. For simplicity, we stick to the
    // linear schedule.

    // For the Ising model, we choose the initial Hamiltonian to be just the
    // uniform transverse field coupling, and the target Hamiltonian to be
    // just the uniform two-site ZZ coupling.

    // Hᵢ = - h ∑ₖ Xₖ,
    // Hₜ = - j ∑ₖ ZₖZₖ₊₁

    // Thus the ground state of Hᵢ is simply the |+〉 product state. The ground
    // state of Hₜ in this case is actually also easy to prepare, but suffices
    // to demonstrate the procedure of adiabatic state preparation.

    // We provide two equivalent solutions to this problem.

    // The first solution manually varies the coefficients on the Ising model
    // GeneratorSystem constructed previously to replicate the interpolated
    // Hamiltonian, which is then packaged as an `EvolutionSchedule` type.
    // This is then fed into the time-dependent simulation algorithm which is
    // of type `TimeDependentSimulationAlgorithm`, and acts on input qubits.

    // The second solution could be more convenient in certain cases. We
    // construct the start Hamiltonian Hᵢ and the target Hamiltonian Hₜ as
    // separate `EvolutionGenerator` types. Together with a choice of
    // `TimeDependentSimulationAlgorithm`, these are then arguments of the
    // library function `AdiabaticEvolution' which automatically interpolates
    // between these Hamiltonians and constructs the `EvolutionSchedule` type
    // and implements time-dependent simulation on the input qubits.

        /// # Summary
    /// We now choose uniform coupling coefficients, allocate qubits to the
    /// simulation, implement adiabatic state preparation, and then return
    /// the results of spin measurement on each site.
    ///
    /// # Input
    /// ## nSites
    /// Number of qubits that the represented system will act upon.
    /// ## hXInitial
    /// Value of the coefficient `h` at s=0.
    /// ## jFinal
    /// Value of the coefficient `j` at s=1.
    /// ## adiabaticTime
    /// Time over which the schedule parameter is varied.
    /// from 0 to 1.
    /// ## trotterStepSize
    /// Time simulated by each step of simulation algorithm.
    /// ## trotterOrder
    /// Order of Trotter–Suzuki integrator.
    ///
    /// # Output
    /// A `Result[]` storing the outcome of Z basis measurements on each site
    /// of the Ising model.
    operation GeneralIsingAdiabaticAndMeasureBuiltIn (nSites : Int, hXInitial : Double, hZFinal : Double[], jZZFinal : Double[], 
	adiabaticTime : Double, trotterStepSize : Double, trotterOrder : Int, qpeStepSize : Double, nBitsPrecision: Int) : Result[] {
        
        let hXCoupling = UniformHCoupling(hXInitial, _);

        // For antiferromagnetic coupling, choose jFinal to be negative.
        let hZCouplingFinal = HCoupling(hZFinal, _);
        let jZZCoupling = JZZCoupling(nSites, jZZFinal, _);

        using (qubits = Qubit[nSites]) {
            Prepare1DIsingState(qubits);
            (IsingAdiabaticEvolutionBuiltIn(
                nSites, adiabaticTime, trotterStepSize, trotterOrder, hXCoupling, hZCouplingFinal, jZZCoupling
            ))(qubits);
            
            DumpRegister("..\\..\\..\\QubitsDump.txt", qubits);
            return ForEach(MResetZ, qubits);
        }
    }

    /// # Summary
    /// This initializes the qubits in an easy-to-prepare eigenstate of the
    /// initial Hamiltonian.
    ///
    /// # Input
    /// ## qubits
    /// Qubit register encoding the Ising model quantum state.
    operation Prepare1DIsingState(qubits : Qubit[]) : Unit is Adj + Ctl {
        ApplyToEachCA(H, qubits);
    }
    
    // We now define functions for the coefficients

    /// # Summary
    /// A function that outputs uniform single-site coupling coefficients
    /// `hₖ`.
    ///
    /// # Input
    /// ## amplitude
    /// Value of coefficient.
    /// ## idxQubit
    /// Index `k` of the qubit that the represented term will act upon.
    ///
    /// # Output
    /// A function returning coefficients `hₖ` for each site.
    function UniformHCoupling(amplitude : Double, idxQubit : Int) : Double {
        return 1.0 * amplitude;
    }

    function HCoupling(h : Double[], idxQubit : Int) : Double {
        return h[idxQubit];
    }

    function JZZCoupling(numQubits : Int, J : Double[], idx_k : Int) : Double {
        let idxQubits = GetTwoSiteQubitIndicesFromK(numQubits, idx_k);

        return J[idxQubits[0] * numQubits + idxQubits[1]];
    }

    /// `TimeDependentSimulationAlgorithm`, and the time of simulation.
    ///
    /// # Input
    /// ## nSites
    /// Number of qubits that the represented system will act upon.
    /// ## adiabaticTime
    /// Time over which the schedule parameter is varied.
    /// from 0 to 1.
    /// ## trotterStepSize
    /// Time simulated by each step of simulation algorithm.
    /// ## trotterOrder
    /// Order of Trotter–Suzuki integrator.
    /// ## hXCoupling
    /// Function returning coefficients `hₖ` for each site.
    /// ## jCoupling
    /// Function returning coefficients `jₖ` for each two-site interaction.
    ///
    /// # Output
    /// a Unitary operator implementing time-dependent evolution by the
    /// Hamiltonian H(s) when s is varied uniformly between 0 and 1 over Time
    /// `adiabaticTime`.
    function IsingAdiabaticEvolutionBuiltIn (nSites : Int, adiabaticTime : Double, trotterStepSize : Double, trotterOrder : Int, hXCoupling : (Int -> Double), hZCouplingFinal : (Int -> Double), jZZCouplingFinal : (Int -> Double)) : (Qubit[] => Unit is Adj + Ctl) {

        // This is the initial Hamiltonian
        let start = StartEvoGen(nSites, hXCoupling);

        // This is the final Hamiltonian
        let end = EndEvoGen(nSites, hZCouplingFinal, jZZCouplingFinal);

        // We choose the time-dependent Trotter–Suzuki decomposition as
        // our simulation algorithm.
        let timeDependentSimulationAlgorithm = TimeDependentTrotterSimulationAlgorithm(trotterStepSize, trotterOrder);

        // The function InterpolatedEvolution uniformly interpolates between the start and the end Hamiltonians.
        return InterpolatedEvolution(adiabaticTime, start, end, timeDependentSimulationAlgorithm);
    }


    //////////////////////////////////////////////////////////////////////////
    // Time-dependent simulation using more built-in functions ///////////////
    //////////////////////////////////////////////////////////////////////////

    // However, in some cases, we are provided with a description of the
    // both Hamiltonians separately, and would like to avoid the need to
    // manually implement this interpolation. A complete description of
    // a Hamiltonian is an `EvolutionGenerator` type that contains
    // both a `GeneratorSystem` that describes terms, and an `EvolutionSet`
    // that maps each term to time-evolution by that term. This will be our
    // starting point.

    /// # Summary
    /// This specifies the initial and target Hamiltonians as separate
    /// `EvolutionGenerator` types.
    ///
    /// # Input
    /// ## nSites
    /// Number of qubits that the represented system will act upon.
    /// ## hXCoupling
    /// Function returning coefficients `hₖ` for each site.
    /// ## jCoupling
    /// Function returning coefficients `jₖ` for each two-site interaction.
    ///
    /// # Output
    /// A `EvolutionGenerator` representing time evolution by each term of the
    /// initial and target Hamiltonians respectively.
    function StartEvoGen (nSites : Int, hXCoupling : (Int -> Double)) : EvolutionGenerator {
        let XGenSys = OneSiteGeneratorSystem(1, nSites, hXCoupling);
        return EvolutionGenerator(PauliEvolutionSet(), XGenSys);
    }

    function EndEvoGen (nSites : Int, hZCoupling : (Int -> Double), jZZCoupling : (Int -> Double)) : EvolutionGenerator {
        let ZGenSys = OneSiteGeneratorSystem(3, nSites, hZCoupling);
        let ZZGenSys = TwoSiteGeneratorSystem(3, nSites, jZZCoupling);
        let combinedGenSys = AddGeneratorSystems(ZGenSys, ZZGenSys);
        return EvolutionGenerator(PauliEvolutionSet(), combinedGenSys);
    }

    //////////////////////////////////////////////////////////////////////////
    // Introduction //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    // In this sample, we demonstrate use of the generator representation
    // functionality offered by the Q# canon to represent Ising model
    // Hamiltonians.

    // Later, we will extend these techniques to represent the
    // 1D Heisenberg XXZ model.

    // We will begin by constructing a representation of the 1D transverse
    // Ising model Hamiltonian,
    //     H = - ( J₀ Z₀ Z₁ + J₁ Z₁ Z₂ + … ) - (h₀ X₀ + h₁ X₁ + …),
    // where {Jᵢ} are nearest-neighbor couplings, and where hₓ is a
    // transverse field.

    // Since this Hamiltonian is naturally expressed in the Pauli basis,
    // we will use the PauliEvolutionSet() function to obtain a simulatable
    // basis to use in representing H. Thus, we begin by defining our
    // indices with respect to the Pauli basis. In doing so, we will define
    // helper functions to return single-site and two-site generator indices.

    //////////////////////////////////////////////////////////////////////////
    // 1D Ising model ////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    /// # Summary
    /// Returns a generator index that is supported on a single site.
    ///
    /// # Input
    /// ## idxPauli
    /// Index of the Pauli operator to be represented, where `1` denotes
    /// `PauliY` and `3` denotes `PauliZ`.
    /// ## idxQubit
    /// Index `k` of the qubit that the represented term will act upon.
    /// ## hCoupling
    /// Function returning coefficients `hₖ` for each site. E.g.:
    /// should return the coefficient for the index at `idxQubit = 3`.
    ///
    /// # Output
    /// A `GeneratorIndex` representing the term - hₖ {Xₖ, Yₖ, Zₖ}, where hₖ is the
    /// function `hCoupling` evaluated at the site index `k`, and where
    /// {Xₖ, Yₖ, Zₖ}, selected by idxPauli, is the Pauli operator acting at the
    /// site index `k`.
    function OneSiteGeneratorIndex (idxPauli : Int, idxQubit : Int, hCoupling : (Int -> Double)) : GeneratorIndex {
        let coeff = 1.0 * hCoupling(idxQubit);
        let idxPauliString = [idxPauli];
        let idxQubits = [idxQubit];
        return GeneratorIndex((idxPauliString, [coeff]), idxQubits);
    }


    /// # Summary
    /// Returns a generator system for a sum of generator indices each
    /// supported on a single site.
    ///
    /// # Input
    /// ## idxPauli
    /// Index of the Pauli operator to be represented.
    /// ## nSites
    /// Number of qubits that the represented system will act upon.
    /// ## hCoupling
    /// Function returning coefficients `hₖ` for each site.
    ///
    /// # Output
    /// A `GeneratorSystem` representing the sum - Σₖ hₖ {Xₖ, Yₖ, Zₖ}.
    function OneSiteGeneratorSystem (idxPauli : Int, nSites : Int, hCoupling : (Int -> Double)) : GeneratorSystem {
        return GeneratorSystem(nSites, OneSiteGeneratorIndex(idxPauli, _, hCoupling));
    }

    /// # Summary
    /// Returns a generator index that is supported on two sites.
    ///
    /// # Input
    /// ## idxPauli
    /// Index of the Pauli operator to be represented, where `1` denotes
    /// `PauliY` and `3` denotes `PauliZ`.
    /// ## nSites
    /// Number of qubits that the represented system will act upon.
    /// ## idx_k
    /// Index `k` is the composite index `k = nSites . i + j`
    /// such that `J_{i,j}` can be calculated from `k`. 
    /// ## jCoupling
    /// Function returning coefficients `J_{i,j}` for each two-site interaction.
    /// E.g.: `jCoupling(3, 7)` should return the coefficient for the coupling
    /// between `idxQubit_i = 3` and `idxQubit_j = 7`.
    ///
    /// # Output
    /// A `GeneratorIndex` representing the term - J_{i,j} {X_iX_j, Y_iY_j, Z_iZ_j},
    /// where J_{ij} is the function `jCoupling` evaluated at the site index `i, j`,
    /// and where {X_iX_j, Y_iY_j, Z_iZ_j}, selected by idxPauli, is the Pauli
    /// operator acting between those sites.
    function TwoSiteGeneratorIndex (idxPauli : Int, nSites : Int, idx_k : Int, jCoupling : (Int -> Double)) : GeneratorIndex {
        let idxQubits = GetTwoSiteQubitIndicesFromK(nSites, idx_k);
        let coeff = -1.0 * jCoupling(idx_k);
        let idxPauliString = [idxPauli, idxPauli];
        let generatorIndex = GeneratorIndex((idxPauliString, [coeff]), idxQubits);

        if (idxQubits[0] >= nSites or idxQubits[1] >= nSites) {
            fail "Qubit index must be smaller than number of sites.";
        }

        return generatorIndex;
    }


    /// # Summary
    /// Returns a generator system for a sum of generator indices each
    /// supported on two neighboring sites.
    ///
    /// # Input
    /// ## idxPauli
    /// Index of the Pauli operator to be represented.
    /// ## nSites
    /// Number of qubits that the represented system will act upon.
    /// ## idxQubit
    /// Index `k` of the qubit that the represented term will act upon.
    /// ## jCoupling
    /// Function returning coefficients `Jₖ` for each two-site interaction.
    ///
    /// # Output
    /// A `GeneratorSystem` representing the sum - Σₖ Jₖ{XₖXₖ₊₁, YₖYₖ₊₁, ZₖZₖ₊₁}.
    function TwoSiteGeneratorSystem(idxPauli : Int, nSites : Int, jCoupling : (Int -> Double)) : GeneratorSystem {
        return GeneratorSystem(nSites * (nSites - 1) / 2, TwoSiteGeneratorIndex(idxPauli, nSites, _, jCoupling));
    }

    function GetTwoSiteQubitIndicesFromK(nSites : Int, idx_k : Int) : Int[] {
        mutable kCounter = idx_k;
        mutable numColumnsInRow = nSites - 1;
        mutable idxQubit_i = -1;
        mutable idxQubit_j = 0;
        while (kCounter >= 0) {
            set idxQubit_i = idxQubit_i + 1;
            set idxQubit_j = idxQubit_i + kCounter + 1;
            set kCounter = kCounter - numColumnsInRow;
            set numColumnsInRow = numColumnsInRow - 1;
        }
        
        if (idxQubit_i >= nSites or idxQubit_i >= nSites) {
            fail "Qubit index must be smaller than number of sites.";
        }

        return [idxQubit_i, idxQubit_j];
	}
}


