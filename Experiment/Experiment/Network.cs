using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace Experiment
{
    public class Network<T>
    {
        public static Network<T> Create(int[] topology, Func<int, T[], double> inputTransform, double bias0 = 0, double weight0 = 1 )
        {
            int maxDepth = topology.Length;
            if (maxDepth < 2) throw new ArgumentException("Depth must be at least two levels");

            Dictionary<int, HashSet<Neuron>> neurons = new Dictionary<int, HashSet<Neuron>>();
            HashSet<Input> inputs = new HashSet<Input>();

            for (int l = 0; l < maxDepth; l++)
            {
                HashSet<Neuron> nls;
                if (!neurons.TryGetValue(l, out nls))
                {
                    nls = new HashSet<Neuron>();
                    neurons.Add(l, nls);
                }
                if (l == 0)
                {
                    for (int i = 0; i < topology[l]; i++)
                    {
                        Input input = new Input(string.Format("{0}:{1}", l + 1, i + 1))
                        {
                            Bias = bias0
                        };
                        nls.Add(input);
                        inputs.Add(input);
                    }
                }
                else if (l == maxDepth - 1)
                {
                    for (int i = 0; i < topology[l]; i++)
                    {
                        Output output = new Output(string.Format("{0}:{1}", l+1, i+1))
                        {
                            Bias = bias0
                        };
                        nls.Add(output);
                        HashSet<Neuron> prev = neurons[l - 1];
                        foreach (var pn in prev)
                        {
                            Synapse s = new Synapse(pn, output);
                            s.Weight = weight0;
                            pn.Axons.Add(s);
                            output.Dendrites.Add(s);
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < topology[l]; i++)
                    {
                        Neuron neuron = new Neuron(string.Format("{0}:{1}", l+1, i+1))
                        {
                            Bias = bias0
                        };
                        nls.Add(neuron);
                        HashSet<Neuron> prev = neurons[l - 1];
                        foreach (var pn in prev)
                        {
                            Synapse s = new Synapse(pn, neuron);
                            s.Weight = weight0;
                            pn.Axons.Add(s);
                            neuron.Dendrites.Add(s);
                        }
                    }
                }
            }

            return new Network<T>(inputTransform, inputs.ToArray()) { _maxDepth = maxDepth };
        }

        /// <summary>
        /// Create a network from a string which describes the neurons, connections and weights
        /// </summary>
        /// <param name="networkMap">of the form L:I:B [W] L:I:B ... where L = neuron level, starting with 1, 
        /// I = neuron level index, starting with 1, 
        /// B = value bias for the neuron and 
        /// W = connection weight</param>
        /// <param name="maxDepth"></param>
        /// <param name="inputTransform"></param>
        /// <returns></returns>
        public static Network<T> Create(string networkMap, int maxDepth, Func<int, T[], double> inputTransform)
        {
            Dictionary<string, Neuron> neurons = new Dictionary<string, Neuron>();

            string pat = @"(?=((?<connection>(?<source>\d:\d(:(?<sourceBias>[-\d\.]*))*)\s*\[(?<weight>[-\d\.]*)\]\s*(?<dest>\d:\d(:(?<destBias>[-\d\.]*))*))))";
            Regex r = new Regex(pat);
            Match m = r.Match(networkMap);
            while (m.Success)
            {
                string source = m.Groups["source"].Value;
                string dest = m.Groups["dest"].Value;
                string weight = m.Groups["weight"].Value;
                double sourceBias = 0, destBias = 0;
                if (!string.IsNullOrEmpty(m.Groups["sourceBias"].Value))
                {
                    sourceBias = double.Parse(m.Groups["sourceBias"].Value);
                    source = source.Replace(":" + m.Groups["sourceBias"].Value, "");
                }
                if (!string.IsNullOrEmpty(m.Groups["destBias"].Value))
                {
                    destBias = double.Parse(m.Groups["destBias"].Value);
                    dest = dest.Replace(":" + m.Groups["destBias"].Value, "");
                }

                Neuron sourceN = null, destN = null;

                if (!neurons.TryGetValue(source, out sourceN))
                {
                    if (source.StartsWith("1"))
                    {
                        sourceN = new Input(source);
                    }
                    else
                    {
                        sourceN = new Neuron(source);
                    }
                    sourceN.Bias = sourceBias;
                    neurons.Add(source, sourceN);
                }

                if (!neurons.TryGetValue(dest, out destN))
                {
                    if (dest.StartsWith(maxDepth + ":"))
                    {
                        destN = new Output(dest);
                    }
                    else
                    {
                        destN = new Neuron(dest);
                    }

                    destN.Bias = destBias;
                    neurons.Add(dest, destN);
                }

                var s = new Synapse(sourceN, destN) {Weight = double.Parse(weight)};
                sourceN.Axons.Add(s);
                destN.Dendrites.Add(s);

                m = m.NextMatch();
            }

            return new Network<T>(inputTransform, neurons.Values.Where(n => n is Input).ToArray()) { _maxDepth = maxDepth };
        }

        private Neuron[] _inputs;
        private Neuron[] _outputs;
        private Func<int, T[], double> _inputTransform;
        private int _maxDepth;
        public Network(Func<int, T[], double> inputTransform, params Neuron[] inputs)
        {
            _inputs = inputs;
            var outputs = new List<Neuron>();
            Traverse((n, level) => true, (n, level) =>  n is Output, outputs.Add);
            _outputs = outputs.ToArray();
            _inputTransform = inputTransform;
        }

        public Neuron[] Inputs { get { return _inputs; } }
        public double[] Outputs { get; set; }

        public double[] Update(T[] inputSignal)
        {
            SetInputs(inputSignal);
            Propagate((n) =>
            {
                n.Value = 0;
                foreach (var input in n.Dendrites)
                {
                    input.Receiver.Value += input.Signaler.Activation() * input.Weight;
                }
            });
            return ReturnOutputs();
        }

        private double[] ReturnOutputs()
        {
            double[] outputs = new double[_outputs.Length];
            for (int i = 0; i < _outputs.Length; i++)
            {
                Output n = _outputs[i] as Output;
                n.Value = 0;
                foreach (var input in n.Dendrites)
                {
                    input.Receiver.Value += input.Signaler.Activation() * input.Weight;
                }

                outputs[i] = _outputs[i].Activation();
            }
            return outputs;
        }

        private void Propagate(Action<Neuron> action,  bool allowMultipleVisits = true)
        {
            Propagate(action, (n, level) => !(n is Input), allowMultipleVisits);
        }

        private void Propagate(Action<Neuron> action, Func<Neuron, int, bool> actionSelector, bool allowMultipleVisits = true)
        {
            Traverse((n, level) => level <= _maxDepth, actionSelector, action, allowMultipleVisits);
        }

        private void BackPropagate(Action<Neuron> action, bool allowMultipleVisits = true)
        {
            Traverse((n, level) => level >= 1, (n, level) => !(n is Output), action, allowMultipleVisits, false);
        }

        private void Traverse(Func<Neuron, int, bool> traversalSelector, Func<Neuron, int, bool> actionSelector, Action<Neuron> action, bool allowMultipleVisits = false, bool forward = true)
        {
            var visited = new HashSet<Neuron>();
            int level = forward ? 1 : _maxDepth;
            Neuron[] start = forward ? _inputs : _outputs;
            foreach (var n in start)
            {
                _traverseNode = n;
                if (traversalSelector(n, level))
                    TraverseRecurse(n, traversalSelector, actionSelector, action, allowMultipleVisits, forward, ref level, ref visited);
            }
        }

        private Neuron _traverseNode = null;
        private void TraverseRecurse(Neuron n, Func<Neuron, int, bool> traversalSelector, Func<Neuron, int, bool> actionSelector, Action<Neuron> action, bool allowMultipleVisits, bool forward, ref int level, ref HashSet<Neuron> visited)
        {
            if (allowMultipleVisits || !visited.Contains(n))
            {
                if (actionSelector(n, level))
                {
                    action(n);
                }
                visited.Add(n);
            }
            
            if (forward)
            {
                foreach (var nn in n.Axons)
                {
                    level++;
                    if (traversalSelector(nn.Receiver, level))
                    {
                        _traverseNode = n;
                        TraverseRecurse(nn.Receiver, traversalSelector, actionSelector, action, allowMultipleVisits, forward, ref level, ref visited);
                        _traverseNode = n;
                    }
                    level--;
                }
            }
            else
            {
                foreach (var nn in n.Dendrites)
                {
                    level--;
                    if (traversalSelector(nn.Signaler, level))
                    {
                        _traverseNode = n;
                        TraverseRecurse(nn.Signaler, traversalSelector, actionSelector, action, allowMultipleVisits, forward, ref level, ref visited);
                        _traverseNode = n;
                    }
                    level++;
                }
            }
        }

        protected void SetInputs(T[] inputSignal)
        {
            _inputIndex = 0;
            Traverse((n, level) => n is Input, (n, level) => true, (n) => SetInputs(n, inputSignal), false);
        }

        private int _inputIndex = 0;
        protected virtual void SetInputs(Neuron n, T[] inputSignal)
        {
            n.Value = _inputTransform(_inputIndex, inputSignal);
            _inputIndex++;
        }

        public virtual void Train(T[][] inputSignalsSet,
            double[][] targetsSet,
            bool reTrain = true,
            double eta = 0.05,
            double alpha = 0.6,
            double w0 = 0.2,
            double epochMax = 100000,
            double errorMin = 0.00002,
            Func<double, double, double> errorFunction = null)
        {
            if (!inputSignalsSet.Length.Equals(targetsSet.Length)) throw new ArgumentException("Input signal set and target output set must be the same length.");

            if (errorFunction == null)
            {
                errorFunction = (output, target) =>
                {
                    double delta = target - output;
                    double error = 0.5 * delta * delta;
                    return error;
                };
            }

            int inputSignalCount = inputSignalsSet.Length;
            
            if (reTrain)
                InitializeWeightsAndDeltas(w0);

            Dictionary<object, double> deltas = new Dictionary<object, double>();
            double perror = 0;
            double derror = 0;
            double pderror = 0;
            double d2error = 0;

            for (int epoch = 0; epoch < epochMax; epoch++)
            {
                T[][] inputSignalsSetR;
                double[][] targetsSetR;
                double error = 0;
                Randomize(inputSignalsSet, targetsSet, out inputSignalsSetR, out targetsSetR);
                

                for (int p = 0; p < inputSignalCount; p++)
                {
                    T[] inputSignals = inputSignalsSetR[p];
                    double[] targets = targetsSetR[p];

                    double[] outputs = Update(inputSignals);

                    for (int o = 0; o < outputs.Length; o++)
                    {
                        ((Output)_outputs[o]).Target = targets[o];
                        double dErr = errorFunction(outputs[o], targets[o]);
                        error += dErr;
                    }

                    for (int i = 0; i < _inputs.Length; i++)
                        _inputs[i].Delta(); // compute and store delta values in entire network

                    HashSet<object> props = new HashSet<object>();
                    Propagate((n) =>
                    {
                        double delta = n.LastDelta;
                        if (!props.Contains(n)) // don't adjust bias more than once during network traversal
                        {
                            double dWN;
                            if (!deltas.TryGetValue(n, out dWN))
                                deltas.Add(n, 0);
                            dWN = eta * delta + alpha * dWN;
                            deltas[n] = dWN;
                            n.Bias += dWN;
                            props.Add(n);
                        }
                        foreach (var nD in n.Dendrites)
                        {
                            if (!props.Contains(nD)) // don't adjust weights more than once during network traversal
                            {
                                double dWND;
                                if (!deltas.TryGetValue(nD, out dWND))
                                    deltas.Add(nD, 0);
                                dWND = eta * nD.Signaler.LastActivation * delta + alpha * dWND;
                                nD.Weight += dWND;
                                deltas[nD] = dWND;
                                props.Add(nD);
                            }
                        }
                    }, true);
                }
                if (error <= errorMin) return;
                if (epoch % 10000 == 0) Console.WriteLine("epoch: {0}; error: {1}; derror: {2}, d2error: {3}", epoch, error, derror, d2error);

                if (epoch > 0) derror = error - perror;
                if (epoch > 1)
                {
                    d2error = derror - pderror;
                }
                pderror = derror;
                perror = error;
                error = 0;
            }
        }

        private void InitializeWeightsAndDeltas(double w0)
        {
            Random r = new Random();
            Propagate((n) =>
            {
                n.Bias = 0.0;
                if (!(n is Input)) n.Value = 0;
                foreach (var dend in n.Dendrites)
                    dend.Weight = w0 * (2 * (r.NextDouble() - 0.5));
            }, (n, l) => true, false);
        }

        private void Randomize(T[][] signals, double[][] targets, out T[][] signalsR, out double[][] targetsR)
        {
            List<T[]> randS = new List<T[]>();
            List<T[]> sourceS = new List<T[]>(signals);
            List<double[]> randT = new List<double[]>();
            List<double[]> sourceT = new List<double[]>(targets);
            Random r = new Random();
            while (randS.Count < signals.Length)
            {
                int i = r.Next(0, sourceS.Count - 1);
                randS.Add(sourceS[i]);
                sourceS.RemoveAt(i);
                randT.Add(sourceT[i]);
                sourceT.RemoveAt(i);
            }
            signalsR = randS.ToArray();
            targetsR = randT.ToArray();
        }

        private Random _r = new Random();
        private double Gaussian(double mean, double stdDev)
        {
            return mean 
                + stdDev 
                * (Math.Sqrt(-2.0 * Math.Log(_r.NextDouble()))) 
                * (Math.Sin(2.0 * Math.PI * _r.NextDouble()));
        }

        /// <summary>
        /// Uses hybrid Evolutionary Algorithm for determining optimized network topology and parameters
        /// </summary>
        /// <param name="signals">training input signals</param>
        /// <param name="targets">training expected outputs</param>
        /// <param name="generations">maximum number of evolutions to execute</param>
        /// <param name="S">Maximum consecutive generations allowable for GL >= 0</param>
        /// <param name="Ps">Structural mutation probability - determines liklihood of changing network topology</param>
        /// <param name="Pw">Weight mutation probability - determines whether network weights will be mutated with random Gaussian noise of step size SS</param>
        /// <param name="SS">Step size of weight mutation</param>
        /// <param name="Pgs">Global structural mutation probability - applies to entire population</param>
        /// <param name="Pgw">Global weight mutation probability - applies to entire population</param>
        /// <param name="SSg">Global step size of weight mutation - applies to entire population</param>
        /// <param name="Pls">Local structural mutation probability - applies to individual</param>
        /// <param name="Plw">Local weight mutation probability - applies to individual</param>
        /// <param name="SSl">Local step size of weight mutation - applies to individual</param>
        /// <param name="tPgs">Temporary global structural mutation probability - determines decay rate of Pgs</param>
        /// <param name="tPgw">Temporary global weight mutation probability - determines decay rate of Pgw</param>
        /// <param name="tSSg">Temporary global step size of weight mutation probability - determines decay rate of SSg</param>
        /// <param name="bPls">Base local structural mutation probability - determines max value of Pls</param>
        /// <param name="bPlw">Base local weight mutation probability - determines max value of Plw</param>
        /// <param name="bSSl">Base local step size of weight mutation probability - determines max value of SSl</param>
        public virtual Network<T> Evolve(
            T[][] signals, 
            double[][] targets,
            int maxHiddenNeurons = 50,
            int populationSize = 100,
            int generations = 10000,
            int S = 500,
            double bPls = 0.2,
            double bPlw = 0.2,
            double bSSl = 2,
            double a1 = 1,
            double a2 = 0.1,
            double a3 = 1)
        {
            double tPgs = bPls;
            double tPgw = bPlw;
            double tSSg = bSSl;
            int sCount = S;
            Func<double> s = () => sCount;
            Func<double> Pgs = () => tPgs * s()/S;    // (1)
            Func<double> Pgw = () => tPgw * s()/S;    // (2)
            Func< double> SSg = () => tSSg * s()/S;    // (3)
            //Func<Func<double>, Func<double>, double> GL = (eb, eo) => ((eb() + 1)/(eo() + 1)) - 1;

            #region algorithm
            /*
             * FORMULAE
             * 
             * (1)  Pgs(t) = tmpPgs ∗ s(t) / S
             * (2)  Pgw(t) = tmpPgw ∗ s(t) / S
             * (3)  SSg(t) = tmpSSg ∗ s(t) / S
             * (4)  Pls(n) = basePls ∗ [1 − Fbest/F(n)]
             * (5)  Plw(n) = basePlw ∗ [1 − Fbest/F(n)]
             * (6)  SSl(n) = baseSSl ∗ [1 − Fbest/F(n)]
             * (7)  Ps(n,t) = Pgs(t) + Pls(n)
             * (8)  Pw(n,t) = Pgw(t) + Plw(n)
             * (9)  SS(n,t) = SSg(t) + SSl(n)
             * (10) wij = w'ij + N(0,(SSl + SSg))
             * 
             * ALGORITHM
             * 1)   Generate an initial population of μ ANNs with
                    random hidden nodes and random connections. Assign the
                    initial weights including the network bias with uniform random
                    values between −1 and 1. Initialize S, basePls, basePlw,
                    baseSSl, tmpPgs, tmpPgw, and tmpSSg.
             * 2)   Calculate the local structural mutation probability
                    (Pls), local weight mutation probability (Plw), and local step
                    size of weight perturbation (SSl) of each individual network
                    in the population according to (4)–(6), respectively. Pls, Plw,
                    and SSl are used to determine the severity of the mutation in
                    a given individual in the population.
             * 3)   Calculate Pgs(t), Pgw(t), and SSg(t) at generation t
                    according to (1)–(3), respectively. Pgs(t), Pgw(t), and SSg(t)
                    are used to guide the entire population’s exploration and
                    exploitation of the search space. Note that this step is executed
                    at intervals of generations (10 generations are used in the
                    current implementation) instead of every generation.
             * 4)   Perform the structural mutation (adding or deleting
                    hidden nodes and network connections) and weight mutation
                    (Gaussian perturbation of weights) on the node vector and
                    connection weight matrix (Fig. 2) for each individual in the
                    parent population. Hidden nodes of an ANN are mutated
                    by bit flipping in the node vector. Connection deletion is
                    performed by defining the non-zero value in the connection
                    matrix as zero. Conversely, connection addition is achieved
                    by initiating a zero value in the connection matrix with
                    uniform random values between −1 and 1. The probability of
                    structural mutation is based on (7). Next, the probability that
                    a non-zero connection weight of a network will be mutated is
                    based on (8). If successful, this particular connection weight
                    is mutated by adding a Gaussian perturbation with mean 0
                    and a step size as described in (9), which is defined as
                    follows:
                                wij = w'ij + N(0,(SSl + SSg))                       (10)
                    where wij and w
                    ij denote the parent and offspring connection
                    weight values in the connection matrix. N(0,σ) is the
                    Gaussian perturbation with mean 0 and standard deviation σ.
                    In this case, σ = SSl + SSg.
             * 5)   Evaluate each offspring created after the mutation
                    according to the classification error in the validation set. If
                    the best network found in the current generation is better than
                    the reserved best network (ANNbest) up to this generation,
                    replace ANNbest with this network and reset the stopping
                    counter to its initial value, S. Under the same condition, update
                    the temporary global structural mutation probability (tmpPgs),
                    temporary global weight mutation probability (tmpPgw), and
                    temporary global step size of weight perturbation (tmpSSg) to
                    Pgs(t), Pg(t)w, and SSg(t), respectively. If there is a tie in the
                    classification error at the validation set when finding ANNbest,
                    the individual with the lowest classification error at the training
                    set will be the new ANNbest. If a tie still exists, the individual
                    with the fewest connections and input nodes will be the new
                    ANNbest. Decrease the stopping counter by one and do not
                    update tmpPgs, tmpPgw, and tmpSSg if the classification error
                    at the validation set of the best network found in the current
                    generation is equal to or greater than ANNbest.
             * 6)   Stop the evolutionary process and go to Step 8 if
                    the stopping counter is decreased to zero or the maximum
                    number of generations has been reached. Otherwise, proceed
                    to the next step.
             * 7)   Sort the offspring according to their fitness values
                    and use the rank-based selection to choose μ–2 offspring to
                    become parents for the next generation. Elitism is also used
                    in HEANN, where one fittest parent and one fittest offspring
                    from the current generation will be retained for the next
                    generation. Hence, μ ANNs are selected as parents for the
                    next generation. Steps 2–7 are repeated until the termination
                    criterion for stopping the evolutionary process is satisfied.
             * 8)   The best network in terms of the classification error
                    in the validation set (ANNbest) is the final ANN for the given
                    problem. 
            */
            #endregion

            List<Fitness> networks = CreateNetworks(
                signals[0].Length, 
                targets[0].Length, 
                maxHiddenNeurons,
                populationSize);

            Fitness best = null;
            double pgs = 0, pgw = 0, ssg = 0;
            while (generations > 0 && sCount > 0)
            {
                if (generations % 10 == 0 || best == null)
                {
                    // update global mutation rates every 10 generations
                    pgs = Pgs();
                    pgw = Pgw();
                    ssg = SSg();
                }

                if (best != null)
                {
                    // mutate network, keeping top two parents for next gen
                    Mutate(networks, 2, pgs, pgw, ssg);
                }

                // run the networks
                for (int p = 0; p < targets.Length; p++)
                {
                    double[] target = targets[p];
                    T[] signal = signals[p];
                    for (int n = 0; n < populationSize; n++)
                    {
                        networks[n].Network.Update(signal);
                        networks[n].Increment(target);
                    }
                }

                // order by best fit
                networks.Sort((f1, f2) =>
                {
                    int val = -f1.F.CompareTo(f2.F);
                    if (val == 0)
                    {
                        // fitness function values match, so prefer lesser complexity
                        val = -f1.Complexity.CompareTo(f2.Complexity);
                    }
                    return val;
                });

                if (best == null || networks[0].F < best.F)
                {
                    // new best fit in this generation was found
                    best = networks[0];
                    // reset local min escape counter
                    sCount = S;
                    // update global mutation rates
                    tPgs = Pgs();
                    tPgw = Pgw();
                    tSSg = SSg();
                }

                sCount--;
                generations--;
            }

            return best.Network;
        }

        private void Mutate(List<Fitness> networks, int keep, 
            double Pgs, 
            double Pgw, 
            double SSg)
        {
            double Fbest = networks[0].F;
            Func<Fitness, double> Ps = (n) => Pgs + n.Pls(Fbest); // (7)
            Func<Fitness, double> Pw = (n) => Pgw + n.Plw(Fbest); // (8)
            Func<Fitness, double> SS = (n) => SSg + n.SSl(Fbest); // (9)
        }

        private List<Fitness> CreateNetworks(int inputCount, int outputCount, int maxHiddenNeurons, int populationSize)
        {
            throw new NotImplementedException();
        }

        class Fitness
        {
            public Fitness(Network<T> network, double a1, double a2,
                double basePls, double basePlw, double baseSSl)
            {
                Network = network;
                A1 = a1;
                A2 = a2;
                BasePls = basePls;
                BasePlw = basePlw;
                BaseSSl = baseSSl;
                Reset();
            }

            public Network<T> Network { get; private set; }

            public double TrainingError
            {
                get { return (100.0/(_patternCount*Network.Outputs.Length))*_sumSqr; }
            }
            public double Complexity { get; private set; }
            public double A1 { get; private set; }
            public double A2 { get; private set; } 
            public double BasePls { get; private set; }
            public double BasePlw { get; private set; }
            public double BaseSSl { get; private set; }
  

            public double Pls(double Fbest)
            {
                return BasePls*(1 - Fbest/F); // (4)
            }
            public double Plw(double Fbest)
            {
                return BasePlw * (1 - Fbest / F); // (5)
            }
            public double SSl(double Fbest)
            {
                return BaseSSl * (1 - Fbest / F); // (6)
            }

            private bool _updateF = true;
            private double _f = 0.0;
            public double F
            {
                get
                {
                    if (_updateF)
                    {
                        _f = A1*TrainingError + A2*Complexity;
                        _updateF = false;
                    }
                    return _f;
                }
            }

            private int _patternCount = 0;
            private double _sumSqr = 0.0;
            public void Increment(double[] targets)
            {
                _updateF = true;
                double sumSqrT = 0.0;
                for (int o = 0; o < targets.Length; o++)
                {
                    sumSqrT += Math.Pow(targets[o] - Network.Outputs[o], 2);
                }
                _sumSqr += sumSqrT;
                _patternCount++;
            }

            public void Reset()
            {
                Complexity = Math.Log(Network.ConnectionCount);
                _patternCount = 0;
                _sumSqr = 0.0;
                _updateF = true;
            }
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            HashSet<object> written = new HashSet<object>();
            
            Func<Neuron, int, bool> traverse = (n, level) =>
            {
                return true;
            };
            Func<Neuron, int, bool> actionSelect = (n, level) =>
            {
                return true;
            };
            Action<Neuron> action = (n) =>
            {
                Synapse s = n.Dendrites.Where(d => d.Signaler.Equals(_traverseNode)).FirstOrDefault();
                if (s != null)
                {
                    if (written.Contains(s))
                    {
                        sb.Append(string.Format("{0}\n", s.Signaler.Name));
                    }
                    else
                    {
                        written.Add(s);
                        if (written.Contains(s.Signaler))
                        {
                            sb.Append(string.Format("{0}\t\t[{1}]\t", s.Signaler.Name, s.Weight.ToString("#0.0#####")));
                        }
                        else
                        {
                            sb.Append(string.Format("{0}:{1}\t[{2}]\t", s.Signaler.Name, s.Signaler.Bias.ToString("#0.0#####"), s.Weight.ToString("#0.0#####")));
                            written.Add(s.Signaler);
                        }

                        if (n.Axons.Count == 0)
                        {
                            if (written.Contains(n))
                                sb.Append(string.Format("{0}\n", n.Name));
                            else
                            {
                                sb.Append(string.Format("{0}:{1}\n", n.Name, n.Bias.ToString("#0.0#####")));
                                written.Add(n);
                            }
                        }
                    }
                }
            };
            Traverse(traverse, actionSelect, action, true);
            return sb.ToString();
        }

        public double ConnectionCount { get; set; }
    }
}
