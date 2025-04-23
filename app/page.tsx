'use client'
import Image from "next/image";
import { HybridACOPSO } from "./Implementation";
import { useState, useRef, useEffect } from 'react';
import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  PointElement,
  LineElement,
  Title,
  Tooltip,
  Legend,
  ScatterController
} from 'chart.js';
import { Line, Scatter } from 'react-chartjs-2';


ChartJS.register(
  CategoryScale,
  LinearScale,
  PointElement,
  LineElement,
  Title,
  Tooltip,
  Legend,
  ScatterController
);

export default function Home() {
  interface City {
    x: number;
    y: number;
  }

  interface OptimizationResult {
    route: number[] | null; // Updated to allow null
    distance: number;
    history: number[];
    iterationData: { // Added this field from the actual result object
      iteration: number;
      aco_distance: number | any[] | null;
      pso_distance: number | any[] | null;
      best_distance: number;
    }[];
    cities: number[][];
    routeCoordinates: number[][];
    convergenceData: {
      xValues: number[];
      yValues: number[];
      title: string;
      xLabel: string;
      yLabel: string;
    };
    routePlotData: {
      cityCoordinates: number[][];
      routePath: number[][];
      title: string;
    };
  }

  const [cities, setCities] = useState<City[]>([]);
  const [newCityX, setNewCityX] = useState<string>('');
  const [newCityY, setNewCityY] = useState<string>('');
  const [result, setResult] = useState<OptimizationResult | null>(null);
  const [isOptimizing, setIsOptimizing] = useState<boolean>(false);
  const [parameters, setParameters] = useState({
    n_ants: 20,
    n_particles: 20,
    alpha: 1.0,
    beta: 5.0,
    rho: 0.5,
    max_iterations: 100
  });
  
  const canvasRef = useRef<HTMLCanvasElement>(null);

  const inputCity = () => {
    const x = parseFloat(newCityX);
    const y = parseFloat(newCityY);
    
    if (!isNaN(x) && !isNaN(y)) {
      setCities([...cities, { x, y }]);
      setNewCityX('');
      setNewCityY('');
    }
  };

  const removeCity = (index: number) => {
    setCities(cities.filter((_, i) => i !== index));
  };

  const generateRandomCities = (count: number) => {
    const newCities = Array(count).fill(0).map(() => ({
      x: Math.random() * 100,
      y: Math.random() * 100
    }));
    setCities(newCities);
  };

  const clearCities = () => {
    setCities([]);
    setResult(null);
  };

  const runOptimization = () => {
    if (cities.length < 4) {
      alert("Please add at least 4 cities to run the optimization");
      return;
    }

    setIsOptimizing(true);

    // Convert cities to the format expected by HybridACOPSO
    const cityCoordinates = cities.map(city => [city.x, city.y]);

    // Run optimization in the next tick to allow UI update
    setTimeout(() => {
      try {
        const hybrid = new HybridACOPSO({
          cities: cityCoordinates,
          n_ants: parameters.n_ants,
          n_particles: parameters.n_particles,
          alpha: parameters.alpha,
          beta: parameters.beta,
          rho: parameters.rho,
          max_iterations: parameters.max_iterations
        });

        
        const optimizationResult = hybrid.optimize() as OptimizationResult;
        setResult(optimizationResult);
      } catch (error) {
        console.error("Optimization error:", error);
        alert("An error occurred during optimization");
      } finally {
        setIsOptimizing(false);
      }
    }, 100);
  };

  // Draw cities and route on canvas
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    // Clear canvas
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    // Scale factors to fit within canvas
    const padding = 20;
    const maxX = Math.max(...cities.map(city => city.x), 100);
    const maxY = Math.max(...cities.map(city => city.y), 100);
    const scaleX = (canvas.width - 2 * padding) / maxX;
    const scaleY = (canvas.height - 2 * padding) / maxY;

    // Draw cities
    cities.forEach((city, index) => {
      ctx.beginPath();
      ctx.arc(
        padding + city.x * scaleX, 
        padding + city.y * scaleY, 
        5, 0, 2 * Math.PI
      );
      ctx.fillStyle = 'blue';
      ctx.fill();
      ctx.closePath();
      
      // Add city index label
      ctx.fillStyle = 'black';
      ctx.font = '12px Arial';
      ctx.fillText(
        `${index}`, 
        padding + city.x * scaleX + 8, 
        padding + city.y * scaleY + 4
      );
    });

    // Draw route if result exists and route is not null
    if (result && result.route) {
      ctx.beginPath();
      const firstCity = cities[result.route[0]];
      ctx.moveTo(
        padding + firstCity.x * scaleX,
        padding + firstCity.y * scaleY
      );
      
      // Draw route lines
      for (let i = 1; i < result.route.length; i++) {
        const city = cities[result.route[i]];
        ctx.lineTo(
          padding + city.x * scaleX,
          padding + city.y * scaleY
        );
      }
      
      // Close the loop
      ctx.lineTo(
        padding + firstCity.x * scaleX,
        padding + firstCity.y * scaleY
      );
      
      ctx.strokeStyle = 'red';
      ctx.lineWidth = 2;
      ctx.stroke();
    }
  }, [cities, result]);

  // Prepare convergence chart data
  const convergenceChartData = {
    labels: result?.convergenceData?.xValues || [],
    datasets: [
      {
        label: 'Best Distance',
        data: result?.convergenceData?.yValues || [],
        borderColor: 'rgb(75, 192, 192)',
        backgroundColor: 'rgba(75, 192, 192, 0.5)',
      }
    ]
  };

  const convergenceChartOptions = {
    responsive: true,
    plugins: {
      legend: {
        position: 'top' as const,
      },
      title: {
        display: true,
        text: result?.convergenceData?.title || 'Convergence Chart',
      },
    },
    scales: {
      x: {
        title: {
          display: true,
          text: result?.convergenceData?.xLabel || 'Iteration'
        }
      },
      y: {
        title: {
          display: true,
          text: result?.convergenceData?.yLabel || 'Best Distance'
        }
      }
    }
  };

  return (
    <main className="w-full text-black min-h-screen bg-blue-100 p-8">
      <h1 className="text-3xl font-serif mb-6 text-center">
        Travelling Sales Man Problem using - Hybrid Ant Colony Optimization and Particle Swarm Optimization (HybridACOPSO)
      </h1>
      
      <div className="grid grid-cols-1 md:grid-cols-2 gap-8">
        <div className="bg-gray-50 p-6 rounded-lg shadow">
          <h2 className="text-xl font-bold mb-4">Configuration</h2>
          
          <div className="mb-4">
            <h3 className="font-bold mb-2">Add Cities</h3>
            <div className="flex gap-2 mb-2">
              <input
                type="number"
                value={newCityX}
                onChange={(e) => setNewCityX(e.target.value)}
                placeholder="X coordinate"
                className="border p-2 w-32"
              />
              <input
                type="number"
                value={newCityY}
                onChange={(e) => setNewCityY(e.target.value)}
                placeholder="Y coordinate"
                className="border p-2 w-32"
              />
              <button
                onClick={inputCity}
                className="bg-blue-500 text-white p-2 rounded"
              >
                Add City
              </button>
            </div>
            
            <div className="flex gap-2">
              <button
                onClick={() => generateRandomCities(20)}
                className="bg-green-500 text-white p-2 rounded"
              >
                Generate 20 Random Cities
              </button>
              <button
                onClick={clearCities}
                className="bg-red-500 text-white p-2 rounded"
              >
                Clear All
              </button>
            </div>
          </div>
          
          <div className="mb-4">
            <h3 className="font-bold mb-2">Algorithm Parameters</h3>
            <div className="grid grid-cols-2 gap-2">
              <div>
                <label className="block text-sm">Number of Ants</label>
                <input
                  type="number"
                  value={parameters.n_ants}
                  onChange={(e) => setParameters({...parameters, n_ants: parseInt(e.target.value)})}
                  className="border p-2 w-full"
                />
              </div>
              <div>
                <label className="block text-sm">Number of Particles</label>
                <input
                  type="number"
                  value={parameters.n_particles}
                  onChange={(e) => setParameters({...parameters, n_particles: parseInt(e.target.value)})}
                  className="border p-2 w-full"
                />
              </div>
              <div>
                <label className="block text-sm">Alpha (ACO)</label>
                <input
                  type="number"
                  step="0.1"
                  value={parameters.alpha}
                  onChange={(e) => setParameters({...parameters, alpha: parseFloat(e.target.value)})}
                  className="border p-2 w-full"
                />
              </div>
              <div>
                <label className="block text-sm">Beta (ACO)</label>
                <input
                  type="number"
                  step="0.1"
                  value={parameters.beta}
                  onChange={(e) => setParameters({...parameters, beta: parseFloat(e.target.value)})}
                  className="border p-2 w-full"
                />
              </div>
              <div>
                <label className="block text-sm">Rho (Evaporation)</label>
                <input
                  type="number"
                  step="0.1"
                  value={parameters.rho}
                  onChange={(e) => setParameters({...parameters, rho: parseFloat(e.target.value)})}
                  className="border p-2 w-full"
                />
              </div>
              <div>
                <label className="block text-sm">Max Iterations</label>
                <input
                  type="number"
                  value={parameters.max_iterations}
                  onChange={(e) => setParameters({...parameters, max_iterations: parseInt(e.target.value)})}
                  className="border p-2 w-full"
                />
              </div>
            </div>
          </div>
          
          <div className="mb-4">
            <button
              onClick={runOptimization}
              disabled={isOptimizing || cities.length < 4}
              className={`w-full p-3 rounded font-bold ${
                isOptimizing || cities.length < 4
                  ? 'bg-gray-300 cursor-not-allowed'
                  : 'bg-blue-600 text-white hover:bg-blue-700'
              }`}
            >
              {isOptimizing ? 'Optimizing...' : 'Run Optimization'}
            </button>
          </div>
          
          <div>
            <h3 className="font-bold mb-2">City List ({cities.length})</h3>
            <div className="max-h-60 overflow-y-auto border rounded">
              <table className="w-full">
                <thead className="bg-gray-100">
                  <tr>
                    <th className="p-2 text-left">Index</th>
                    <th className="p-2 text-left">X</th>
                    <th className="p-2 text-left">Y</th>
                    <th className="p-2 text-left">Action</th>
                  </tr>
                </thead>
                <tbody>
                  {cities.map((city, index) => (
                    <tr key={index} className="border-t">
                      <td className="p-2">{index}</td>
                      <td className="p-2">{city.x.toFixed(2)}</td>
                      <td className="p-2">{city.y.toFixed(2)}</td>
                      <td className="p-2">
                        <button
                          onClick={() => removeCity(index)}
                          className="text-red-500 hover:text-red-700"
                        >
                          Remove
                        </button>
                      </td>
                    </tr>
                  ))}
                  {cities.length === 0 && (
                    <tr>
                      <td colSpan={4} className="p-2 text-center text-gray-500">
                        No cities added yet
                      </td>
                    </tr>
                  )}
                </tbody>
              </table>
            </div>
          </div>
        </div>
        
        <div>
          <div className="bg-gray-50 p-6 rounded-lg shadow mb-8">
            <h2 className="text-xl font-bold mb-4">City Map & Optimal Route</h2>
            <canvas 
              ref={canvasRef} 
              width={500} 
              height={400} 
              className="w-full border rounded bg-white"
            />
            {result && (
              <div className="mt-2 text-center font-bold text-green-600">
                Best Route Distance: {result.distance.toFixed(2)}
              </div>
            )}
          </div>
          
          {result && (
            <div className="bg-gray-50 p-6 rounded-lg shadow">
              <h2 className="text-xl font-bold mb-4">Convergence Graph</h2>
              <Line 
                data={convergenceChartData} 
                options={convergenceChartOptions}
              />
            </div>
          )}
        </div>
      </div>
      
      {result && (
        <div className="mt-8 bg-gray-50 p-6 rounded-lg shadow">
          <h2 className="text-xl font-bold mb-4">Optimization Results</h2>
          <div className="mb-4">
            <h3 className="font-bold">Best Route:</h3>
            <div className="p-2 bg-white border rounded overflow-x-auto">
              <code>
                {result.route ? `${result.route.join(' → ')} → ${result.route[0]}` : 'No route found'}
              </code>
            </div>
          </div>
          <div>
            <h3 className="font-bold">Final Distance: {result.distance.toFixed(2)}</h3>
          </div>
        </div>
      )}
    </main>
  );
}