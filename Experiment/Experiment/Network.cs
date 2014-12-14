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
                if (epoch % 10000 == 0) Console.WriteLine("epoch: {0}; error: {1}", epoch, error);
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

        public Neuron[] Inputs { get { return _inputs; } }
        public double[] Outputs { get; set; }

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
    }
}
