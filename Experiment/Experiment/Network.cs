using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace Experiment
{
    [Serializable]
    public class Network<T>
    {
        public static Network<T> Create(int[] topology, Func<int, T[], double> inputTransform, double bias0 = 0, double weight0 = 1)
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
                        Output output = new Output(string.Format("{0}:{1}", l + 1, i + 1))
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
                        Neuron neuron = new Neuron(string.Format("{0}:{1}", l + 1, i + 1))
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

            return new Network<T>(inputTransform, neurons.Values.SelectMany(n => n.ToArray()).ToArray()) { _maxDepth = maxDepth };
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

                var s = new Synapse(sourceN, destN) { Weight = double.Parse(weight) };
                sourceN.Axons.Add(s);
                destN.Dendrites.Add(s);

                m = m.NextMatch();
            }

            return new Network<T>(inputTransform, neurons.Values.ToArray()) { _maxDepth = maxDepth };
        }

        private Neuron[] _inputs;
        private Neuron[] _outputs;
        private Func<int, T[], double> _inputTransform;
        private int _maxDepth;
        public Network(Func<int, T[], double> inputTransform, params Neuron[] neurons)
        {
            Neurons = neurons;
            _inputs = neurons.OfType<Input>().ToArray();
            var outputs = new List<Neuron>();
            _outputs = neurons.OfType<Output>().ToArray();
            _inputTransform = inputTransform;
            Results = new double[_outputs.Length];
        }

        public Neuron[] Inputs { get { return _inputs; } }
        public Neuron[] Outputs { get { return _outputs; } }
        public Neuron[] Neurons { get; private set; }

        public double[] Results { get; private set; }

        public double[] Update(T[] inputSignal, int recursionCount = 1)
        {
            ClearValues();
            SetInputs(inputSignal);
            Propagate((n) =>
            {
                n.Value = 0;
                foreach (var input in n.Dendrites)
                {
                    input.Receiver.Value += input.Signaler.Activation() * input.Weight;
                }
            }, true, recursionCount);
            return ReturnOutputs();
        }

        private void ClearValues()
        {
            foreach (var n in Neurons)
            {
                n.Reset();
            }
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
            Results = outputs;
            return outputs;
        }

        private void Propagate(Action<Neuron> action, bool revisit = true, int recursionCount = 0)
        {
            Propagate(action, (n, level) => !(n is Input), revisit, recursionCount);
        }

        private void Propagate(Action<Neuron> action, Func<Neuron, int, bool> actionSelector, bool revisit = true, int recursionCount = 0)
        {
            Traverse((n, level) => true, actionSelector, action, revisit, recursionCount);
        }

        private void BackPropagate(Action<Neuron> action, bool revisit = true, int recursionCount = 0)
        {
            Traverse((n, level) => level >= 1, (n, level) => !(n is Output), action, revisit, recursionCount, false);
        }

        private void Traverse(Func<Neuron, int, bool> traversalSelector, Func<Neuron, int, bool> actionSelector, Action<Neuron> action, bool revisit, int recursionCount = 0, bool forward = true)
        {
            if (recursionCount < 0) throw new ArgumentException("Recursion count must be 0 or higher");
            var visited = new HashSet<Neuron>();
            int level = forward ? 1 : _maxDepth;
            Neuron[] start = forward ? _inputs : _outputs;
            foreach (var n in start)
            {
                _traverseNode = n;
                if (traversalSelector(n, level))
                {
                    n.RecursionCount++;
                    TraverseRecurse(n, traversalSelector, actionSelector, action, revisit, recursionCount, forward, ref level, ref visited);
                    n.RecursionCount--;
                }
            }
        }

        private Neuron _traverseNode = null;
        private void TraverseRecurse(Neuron n, Func<Neuron, int, bool> traversalSelector, Func<Neuron, int, bool> actionSelector, Action<Neuron> action, bool revisit, int recursionCount, bool forward, ref int level, ref HashSet<Neuron> visited)
        {
            try
            {
                //if (new StackTrace().FrameCount > 25) Debugger.Break();

                bool hasVisited = visited.Contains(n);
                if (revisit || !hasVisited)
                {
                    if (actionSelector(n, level))
                    {
                        action(n);
                    }
                }
                else return;

                if (!hasVisited)
                    visited.Add(n);

                if (recursionCount + 1 >= n.RecursionCount)
                {
                    if (forward)
                    {
                        foreach (var nn in n.Axons)
                        {
                            level++;
                            if (traversalSelector(nn.Receiver, level))
                            {
                                _traverseNode = n;
                                n.RecursionCount++;
                                TraverseRecurse(nn.Receiver, traversalSelector, actionSelector, action, revisit, recursionCount, forward, ref level, ref visited);
                                n.RecursionCount--;
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
                                n.RecursionCount++;
                                TraverseRecurse(nn.Signaler, traversalSelector, actionSelector, action, revisit, recursionCount, forward, ref level, ref visited);
                                n.RecursionCount--;
                                _traverseNode = n;
                            }
                            level++;
                        }
                    }
                }
            }
            finally
            {

            }
        }

        protected void SetInputs(T[] inputSignal)
        {
            _inputIndex = 0;
            Traverse((n, level) => n is Input, (n, level) => true, (n) => SetInputs(n, inputSignal), false, 1);
        }

        private int _inputIndex = 0;
        protected virtual void SetInputs(Neuron n, T[] inputSignal)
        {
            n.Value = _inputTransform(_inputIndex, inputSignal);
            _inputIndex++;
        }

        public virtual TrainingResults AutoTrain(T[][] inputSignalsSet,
            double[][] targetsSet,
            bool reTrain = true,
            double etaBestGuess = 0.000005,
            double w0 = 0.2,
            double epochMax = 100000,
            double errorMin = 0.00002,
            Func<double, double, double> errorFunction = null,
            double converganceRateMin = 0.0,
            int converganceRateEpochs = 500,
            int converganceRollingWindow = 150)
        {
            return Train(inputSignalsSet,
                targetsSet,
                reTrain,
                etaBestGuess,
                0.0,
                w0,
                epochMax,
                errorMin,
                errorFunction,
                converganceRateMin,
                converganceRateEpochs,
                converganceRollingWindow,
                true);
        }

        public virtual TrainingResults ManualTrain(T[][] inputSignalsSet,
            double[][] targetsSet,
            bool reTrain = true,
            double eta = 0.000005,
            double alpha = 0.6,
            double w0 = 0.2,
            double epochMax = 100000,
            double errorMin = 0.00002,
            Func<double, double, double> errorFunction = null,
            double converganceRateMin = 0.0,
            int converganceRateEpochs = 500,
            int converganceRollingWindow = 150)
        {
            return Train(inputSignalsSet,
                targetsSet,
                reTrain,
                eta,
                alpha,
                w0,
                epochMax,
                errorMin,
                errorFunction,
                converganceRateMin,
                converganceRateEpochs,
                converganceRollingWindow,
                false);
        }

        protected virtual TrainingResults Train(T[][] inputSignalsSet,
            double[][] targetsSet,
            bool reTrain,
            double eta,
            double alpha,
            double w0,
            double epochMax,
            double errorMin,
            Func<double, double, double> errorFunction,
            double converganceRateMin,
            int converganceRateEpochs,
            int converganceRollingWindow,
            bool autoAdjustEta)
        {
            Console.WriteLine("eta: {0}, alpha: {1}", eta, alpha);
            if (!inputSignalsSet.Length.Equals(targetsSet.Length)) throw new ArgumentException("Input signal set and target output set must be the same length.");

            Dictionary<object, double> deltas = new Dictionary<object, double>();
            double perror = 0;
            double derror = 0;
            double pderror = 0;
            double d2error = 0;
            double pd2error = 0;
            double d3error = 0;
            int epoch = 0;
            List<double> derrors = new List<double>(new double[] { 0d });
            List<double> d2errors = new List<double>(new double[] { 0d, 0d });
            List<double> d3errors = new List<double>(new double[] { 0d, 0d, 0d });
            int rollingAvCount = 2;

            Action dEta = () =>
            {
                var d2MCount = d2errors.Count();
                var dMCount = derrors.Count();
                var d2M = 0d;

                for (int i = d2MCount - 1; i > d2MCount - 1 - rollingAvCount; i--)
                    d2M += d2errors[i];
                d2M /= rollingAvCount;

                var dM = 0d;
                for (int i = dMCount - 1; i > dMCount - 1 - rollingAvCount; i--)
                    dM += derrors[i];
                dM /= rollingAvCount;

                var ratio = d2M / dM;
                var sign = dM / Math.Abs(dM);
                var sp = eta * (1 + (1 - sign * 0.999 * Math.Tanh(Math.Pow(Math.Abs(d2M), 1d / 8d)))) / 2d;
                eta = sp;
            };

            if (errorFunction == null)
            {
                errorFunction = (output, target) =>
                {
                    double delta = target - output;
                    double err = 0.5 * delta * delta;
                    return err;
                };
            }

            int inputSignalCount = inputSignalsSet.Length;

            if (reTrain)
                InitializeWeightsAndDeltas(w0);

            Console.WriteLine("epoch, error, derror, d2error, d3error, delta, eta");
            double error = 0;
            for (epoch = 0; epoch < epochMax; epoch++)
            {
                T[][] inputSignalsSetR;
                double[][] targetsSetR;

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
                    });
                }

                if (error <= errorMin)
                    return new TrainingResults(epoch,
                        error,
                        errorMin,
                        derrors.Average(),
                        converganceRateMin,
                        eta,
                        alpha,
                        TrainingResults.TrainingReason.Converged);

                if (epoch > 0)
                {
                    derror = error - perror;
                    derrors.Add(derror);
                }
                if (epoch > 1)
                {
                    d2error = derror - pderror;
                    d2errors.Add(d2error);
                }
                if (epoch > 2)
                {
                    d3error = d2error - pd2error;
                    d3errors.Add(d3error);
                }

                if (derrors.Count > converganceRollingWindow)
                {
                    derrors.RemoveAt(0);
                    d2errors.RemoveAt(0);
                    d3errors.RemoveAt(0);
                }

                pd2error = d2error;
                pderror = derror;
                perror = error;

                if (epoch > 0 && epoch % converganceRateEpochs == 0)
                    Console.WriteLine(
                        "{0}, {1}, {2}, {3}, {4}, {5}, {6}",
                        epoch,
                        error,
                        derrors.Average(),
                        d2errors.Average(),
                        d3errors.Average(),
                        d2error / derror,
                        eta);

                if (epoch > rollingAvCount)
                {
                    if (derrors.Average() > converganceRateMin
                        && epoch > converganceRateEpochs
                        && epoch % converganceRollingWindow == 0)
                        return new TrainingResults(epoch,
                            error,
                            errorMin,
                            derrors.Average(),
                            converganceRateMin,
                            eta,
                            alpha,
                            TrainingResults.TrainingReason.DidNotConverge);

                    if (autoAdjustEta)
                        dEta();
                }

                error = 0;
            }

            Console.WriteLine(
                        "{0}, {1}, {2}, {3}, {4}, {5}, {6}",
                        epoch,
                        error,
                        derrors.Average(),
                        d2errors.Average(),
                        d3errors.Average(),
                        d2error / derror,
                        eta);

            return new TrainingResults(epoch,
                            error,
                            errorMin,
                            derrors.Average(),
                            converganceRateMin,
                            eta,
                            alpha,
                            TrainingResults.TrainingReason.DidNotConverge
                            | TrainingResults.TrainingReason.EpochMax);
        }

        public class TrainingResults
        {
            [Flags]
            public enum TrainingReason : int
            {
                Converged = 0x0001,
                DidNotConverge = 0x0010,
                EpochMax = 0x1000
            }

            public TrainingResults(int epochs,
                double finalError,
                double targetError,
                double finalErrorRate,
                double targetErrorRate,
                double eta,
                double alpha,
                TrainingReason reason)
            {
                Epochs = epochs;
                FinalError = finalError;
                TargetError = targetError;
                FinalErrorRate = finalErrorRate;
                TargetErrorRate = targetErrorRate;
                Eta = eta;
                Alpha = alpha;
                Reason = reason;
            }

            public int Epochs { get; private set; }
            public double FinalError { get; private set; }
            public double TargetError { get; private set; }
            public double FinalErrorRate { get; private set; }
            public double TargetErrorRate { get; private set; }
            public TrainingReason Reason { get; private set; }

            public double Eta { get; private set; }
            public double Alpha { get; private set; }
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
            }, (n, l) => true, false, 1);
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

        private static Random _r = new Random();
        private double Gaussian(double mean, double stdDev)
        {
            return mean
                + stdDev
                * (Math.Sqrt(-2.0 * Math.Log(_r.NextDouble())))
                * (Math.Sin(2.0 * Math.PI * _r.NextDouble()));
        }

        public static Network<T> Create(
            Func<int, T[], double> inputTransform,
            T[][] signals,
            double[][] targets,
            int maxHiddenNeurons = 50,
            int populationSize = 100,
            double acceptableFitness = 200,
            int S = 500,
            double bPls = 0.2,
            double bPlw = 0.2,
            double bSSl = 2,
            double a1 = 1,
            double a2 = 0.1,
            double a3 = 1)
        {
            return Create(inputTransform, signals, targets, maxHiddenNeurons, populationSize, int.MaxValue,
                acceptableFitness, S, bPls, bPlw, bSSl, a1, a2, a3);
        }

        /// <summary>
        /// Uses hybrid Evolutionary Algorithm for determining optimized network topology and parameters
        /// </summary>
        /// <param name="signals">training input signals</param>
        /// <param name="targets">training expected outputs</param>
        /// <param name="maxHiddenNeurons">maximum number of interneurons created in evolutionary candidates</param>
        /// <param name="populationSize">number of networks to create per generation</param>
        /// <param name="generations">maximum number of evolutions to execute</param>
        /// <param name="S">Maximum consecutive generations allowable for non-improving evolutions</param>
        /// <param name="bPls">Base local structural mutation probability - determines max value of Pls</param>
        /// <param name="bPlw">Base local weight mutation probability - determines max value of Plw</param>
        /// <param name="bSSl">Base local step size of weight mutation probability - determines max value of SSl</param>
        /// <param name="a1">weight of output fitness in error function</param>
        /// <param name="a2">weight of network complexity in error function</param>
        /// <param name="a3"></param>
        public static Network<T> Create(
            Func<int, T[], double> inputTransform,
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
            return Create(inputTransform, signals, targets, maxHiddenNeurons, populationSize, generations,
                double.MaxValue, S, bPls, bPlw, bSSl, a1, a2, a3);
        }

        private static Network<T> Create(
            Func<int, T[], double> inputTransform,
            T[][] signals,
            double[][] targets,
            int maxHiddenNeurons,
            int populationSize,
            int generations,
            double acceptableError,
            int S,
            double bPls,
            double bPlw,
            double bSSl,
            double a1,
            double a2,
            double a3)
        {
            double tPgs = bPls;
            double tPgw = bPlw;
            double tSSg = bSSl;
            int sCount = S;
            Func<double> s = () => sCount;
            // (1) only modify structure when weight evolution fails, or we're evolve for a fixed set of iterations
            //Func<double> Pgs = () => sCount == 1 || double.MaxValue.Equals(acceptableError) ? tPgs * s() / S : 0;    
            Func<double> Pgs = () => tPgs * s() / S;   // (1)
            Func<double> Pgw = () => tPgw * s() / S;    // (2)
            Func<double> SSg = () => tSSg * s() / S;    // (3)
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
                inputTransform,
                signals[0].Length,
                targets[0].Length,
                maxHiddenNeurons,
                populationSize,
                a1, a2, a3, bPls, bPlw, bSSl);

            double bestF = double.MaxValue;
            Fitness best = null;
            double pgs = 0, pgw = 0, ssg = 0;
            int keep = (int)((double)populationSize * .1);
            while (generations > 0 && sCount > 0)
            {
                if (generations % 10 == 0 || best == null)
                {
                    // update global mutation rates every 10 generations
                    pgs = Pgs();
                    pgw = Pgw();
                    ssg = SSg();
                }

                if (networks[0].F < bestF)
                {
                    int lifespan = best == null ? 0 : best.Generations;
                    // new best fit in this generation was found
                    best = networks[0];
                    bestF = best.F;
                    // reset local min escape counter
                    sCount = S;
                    // update global mutation rates
                    tPgs = Pgs();
                    tPgw = Pgw();
                    tSSg = SSg();
                    Console.WriteLine(string.Format("G: {2}, L: {3}, F: {0}, Complexity: {1}", bestF, best.Complexity, generations, lifespan));
                }
                else if (best != null)
                {
                    networks.Insert(0, best);
                    networks.RemoveAt(networks.Count - 1);
                    best.Generations++;
                }

                // mutate network, keeping top two parents for next gen
                Mutate(networks, keep, pgs, pgw, ssg, acceptableError);

                // run the networks
                ApplySignals(networks, signals, targets);

                // order by best fit
                networks.Sort((f1, f2) =>
                {
                    int val = f1.F.CompareTo(f2.F);
                    if (val == 0)
                    {
                        // fitness function values match, so prefer lesser complexity
                        val = f1.Complexity.CompareTo(f2.Complexity);
                    }
                    return val;
                });

                if (networks[0].F <= acceptableError) break;

                sCount--;
                generations--;
            }

            if (networks[0].F <= acceptableError)
            {
                Network<T> theBest = null;
                double f = double.MaxValue;
                foreach (var n in networks.Where(nn => nn.F <= acceptableError))
                {
                    var results = n.Network.AutoTrain(signals, targets, false);
                    var nF = results.FinalError * n.Complexity;
                    if (nF < f)
                    {
                        theBest = n.Network;
                        f = nF;
                    }
                }
                return theBest;
            }

            return networks[0].Network;
        }

        private static void ApplySignals(List<Fitness> networks, T[][] signals, double[][] targets)
        {
            foreach (var n in networks)
            {
                var result = n.Network.AutoTrain(signals, targets, false, 0.000005, 0.2, 500, 2d, null, 0.0, 50, 25);
                if (result.Reason == TrainingResults.TrainingReason.Converged)
                {
                    result = n.Network.AutoTrain(signals, targets, false, result.Eta, 0.2, 5000, 0.001, null, 0.0, 500, 150);
                }
            }
            for (int p = 0; p < targets.Length; p++)
            {
                double[] target = targets[p];
                T[] signal = signals[p];
                for (int n = 0; n < networks.Count; n++)
                {
                    networks[n].Network.Update(signal, 1);
                    networks[n].Increment(target);
                }
            }
        }

        private static void Mutate(List<Fitness> networks,
            int keep,
            double Pgs,
            double Pgw,
            double SSg,
            double targetF = 10000)
        {
            double Fbest = networks[0].F;
            var nc = networks.Count;

            keep = (int)((float)keep + (float)nc * (targetF / Fbest).Clamp(0d, 1d) / 2f);

            for (int n = 0; n < keep; n++)
            {
                networks[n].Generations++;
                networks[n] = networks[n].Clone(); // clone keepers
            }

            Func<Fitness, double> Ps = (n) => Pgs + n.Pls(Fbest); // (7)
            Func<Fitness, double> Pw = (n) => Pgw + n.Plw(Fbest); // (8)
            Func<Fitness, double> SS = (n) => SSg + n.SSl(Fbest); // (9)

            for (int n = networks.Count - 1; n > keep; n--)
            {
                networks[n] = Mutate(networks[n - keep], Ps, Pw, SS); // clone and mutate keepers
                networks[n].Generations = 0;
            }

            #region mutate keepers only
            //int nn = 1;
            //int nnn = 1;
            //networks[0] = networks[0].Clone(); // don't mutate best fit
            //List<Fitness> keepers = new List<Fitness>(networks.Take(keep));
            //while (nn < networks.Count)
            //{
            //    networks[nn] = Mutate(keepers[nnn], Ps, Pw, SS); // mutate keepers only
            //    nn++;
            //    nnn++;
            //    if (nnn >= keep) nnn = 0;
            //}
            #endregion

        }

        private static Fitness Mutate(Fitness fitness,
            Func<Fitness, double> Ps,
            Func<Fitness, double> Pw,
            Func<Fitness, double> SS)
        {
            return new Fitness(
                fitness.Network.Mutate(Ps(fitness), Pw(fitness), SS(fitness)),
                fitness.A1,
                fitness.A2,
                fitness.A3,
                fitness.BasePls,
                fitness.BasePlw,
                fitness.BaseSSl);
        }

        public Network<T> Clone()
        {
            Dictionary<string, Neuron> neurons = new Dictionary<string, Neuron>();
            foreach (var n in this.Neurons)
            {
                Neuron newN;
                if (n is Dead)
                    newN = new Dead(n.Name);
                else if (n is Input)
                    newN = new Input(n.Name);
                else if (n is Output)
                    newN = new Output(n.Name);
                else
                    newN = new Neuron(n.Name);

                newN.Bias = n.Bias;
                newN.Value = n.Value;

                neurons.Add(newN.Name, newN);
            }

            foreach (var n in this.Neurons)
            {
                Neuron source = neurons[n.Name];
                foreach (var s in n.Axons)
                {
                    Neuron dest = neurons[s.Receiver.Name];
                    var sNew = new Synapse(source, dest) { Weight = s.Weight };
                    source.Axons.Add(sNew);
                    dest.Dendrites.Add(sNew);
                }
            }
            var clone = new Network<T>(this._inputTransform, neurons.Values.ToArray());

            return clone;
        }

        private Network<T> Mutate(double Ps, double Pw, double SS)
        {
            Func<double, double> adjust = (ss) => Gaussian(0, ss);
            _cc = _ccI = _ccO = -1;

            Network<T> child;
            do
            {
                Neuron[] neurons = this.Clone().Neurons;
                // do structure changes first
                for (int n = 0; n < neurons.Length; n++)
                {
                    Neuron source = neurons[n];
                    if (Pw > _r.NextDouble())
                    {
                        source.Bias += adjust(SS);
                    }
                    if (source.GetType().Equals(typeof(Neuron))
                        && source.Axons.Count == 0
                        && source.Dendrites.Count == 0)
                    {
                        source = new Dead(source.Name);
                        neurons[n] = source;
                    }

                    for (int nn = 0; nn < neurons.Length; nn++)
                    {
                        double PsR = _r.NextDouble();
                        double PwR = _r.NextDouble();
                        Neuron target = neurons[nn];
                        if (nn == n || target is Input) continue; // neurons can't connect back to inputs
                        var existing = source.Axons.Where(ss => ss.Receiver.Name.Equals(target.Name)).ToArray();
                        if (existing.Length == 0)
                        {
                            if (Ps > PsR) // check to add a connection
                            {
                                // not connected
                                if (target is Dead)
                                {
                                    target = new Neuron(target.Name) { Bias = target.Bias };
                                }
                                var s = new Synapse(source, target) { Weight = adjust(SS) };
                                source.Axons.Add(s);
                                target.Dendrites.Add(s);
                            }
                        }
                        else
                        {
                            foreach (var s in existing)
                            {
                                if (Ps > PsR) // check to kill the connetion
                                {
                                    if (s.Receiver is Dead)
                                    {
                                        // make it alive
                                        target.Dendrites.Remove(s);
                                        source.Axons.Remove(s);
                                        target = new Neuron(target.Name) { Bias = target.Bias };
                                        var newS = new Synapse(source, target);
                                        source.Axons.Add(newS);
                                        target.Dendrites.Add(newS);
                                    }
                                    else
                                    {
                                        // kill it
                                        target.Dendrites.Remove(s);
                                        source.Axons.Remove(s);
                                    }
                                }
                                if (Pw > PwR) // check to modify its weight
                                {
                                    s.Weight += adjust(SS);
                                }
                                PsR = _r.NextDouble();
                                PwR = _r.NextDouble();
                            }
                        }
                        neurons[nn] = target;
                    }
                }
                child = new Network<T>(this._inputTransform, neurons);
            } while (child.ConnectionCount == 0
                || child.InputConnectionCount == 0
                || child.OutputConnectionCount == 0);


            return child;
        }

        private static List<Fitness> CreateNetworks(
            Func<int, T[], double> inputTransform,
            int inputCount,
            int outputCount,
            int maxHiddenNeurons,
            int populationSize,
            double a1,
            double a2,
            double a3,
            double basePls,
            double basePlw,
            double baseSSl)
        {
            List<Fitness> population = new List<Fitness>();
            for (int n = 0; n < populationSize; n++)
            {
                population.Add(
                    new Fitness(
                            Network<T>.Create(new int[] { inputCount, maxHiddenNeurons, outputCount },
                            inputTransform),
                            a1,
                            a2,
                            a3,
                            basePls,
                            basePlw,
                            baseSSl
                        )
                    );
            }
            return population;
        }
        [Serializable]
        class Fitness
        {
            public Fitness(Network<T> network, double a1, double a2, double a3,
                double basePls, double basePlw, double baseSSl)
            {
                Network = network;
                A1 = a1;
                A2 = a2;
                A3 = a1;
                BasePls = basePls;
                BasePlw = basePlw;
                BaseSSl = baseSSl;
                Reset();
            }

            public Network<T> Network { get; private set; }

            public double TrainingError
            {
                get
                {
                    double val = (100.0 / (_patternCount * Network.Outputs.Length)) * _sumSqr;
                    if (double.NaN.Equals(val))
                    {
                        return double.MaxValue;
                    }
                    return val;
                }
            }
            public double Complexity { get { return Math.Log(Network.ConnectionCount); } }
            public double A1 { get; private set; }
            public double A2 { get; private set; }
            public double A3 { get; private set; }
            public double BasePls { get; private set; }
            public double BasePlw { get; private set; }
            public double BaseSSl { get; private set; }
            public int Generations { get; set; }


            public double Pls(double Fbest)
            {
                return BasePls * (1 - Fbest / F); // (4)
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
                        _f = A1 * TrainingError + A2 * Complexity;
                        _updateF = false;
                    }
                    if (_f < 0) Debugger.Break();
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
                    sumSqrT += Math.Pow(targets[o] - Network.Results[o], 2);
                }
                _sumSqr += sumSqrT;
                _patternCount++;
            }

            public void Reset()
            {
                _patternCount = 0;
                _sumSqr = 0.0;
                _updateF = true;
            }

            public Fitness Clone()
            {
                return new Fitness(this.Network.Clone(),
                    this.A1, this.A2, this.A3, this.BasePls, this.BasePlw, this.BaseSSl) { Generations = this.Generations };
            }

            public override string ToString()
            {
                return string.Format("F: {0}, N:{1}", F, Network);
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
            Traverse(traverse, actionSelect, action, true, 1);
            return sb.ToString();
        }

        public void Save(Stream destination)
        {
            using (StreamWriter sw = new StreamWriter(destination))
            {
                foreach (var n in Neurons)
                {
                    sw.Write(n.ToString());
                    sw.Write(Environment.NewLine);
                    foreach (var s in n.Axons)
                    {
                        sw.Write(s.ToString());
                        sw.Write(Environment.NewLine);
                    }
                }
            }
        }

        int _cc = -1;
        public int ConnectionCount
        {
            get
            {
                if (_cc < 0)
                {
                    _cc = Neurons.Sum(n => n.Axons.Count);
                }
                return _cc;
            }
        }

        int _ccI = -1;
        public int InputConnectionCount
        {
            get
            {
                if (_ccI < 0)
                {
                    _ccI = Inputs.Sum(n => n.Axons.Count);
                }
                return _ccI;
            }
        }

        int _ccO = -1;
        public int OutputConnectionCount
        {
            get
            {
                if (_ccO < 0)
                {
                    _ccO = Outputs.Sum(n => n.Dendrites.Count);
                }
                return _ccO;
            }
        }
    }

    public static class doubleEx
    {
        public static double Clamp(this double value, double min, double max)
        {
            if (value > max) return max;
            if (value < min) return min;
            return value;
        }
    }
}
