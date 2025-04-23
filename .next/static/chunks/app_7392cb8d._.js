(globalThis.TURBOPACK = globalThis.TURBOPACK || []).push([typeof document === "object" ? document.currentScript : undefined, {

"[project]/app/Implementation.js [app-client] (ecmascript)": (function(__turbopack_context__) {

var { g: global, __dirname, k: __turbopack_refresh__, m: module, e: exports } = __turbopack_context__;
{
class HybridACOPSO {
    constructor({ cities, n_ants = 10, n_particles = 10, alpha = 1.0, beta = 2.0, rho = 0.5, w = 0.5, c1 = 1.5, c2 = 1.5, q0 = 0.9, max_iterations = 100 }){
        /**
       * Initialize the hybrid ACO-PSO algorithm.
       * 
       * Parameters:
       * -----------
       * cities : array
       *    Coordinates of cities as [x, y] pairs
       * n_ants : int
       *    Number of ants in the colony
       * n_particles : int
       *    Number of particles in the swarm
       * alpha : float
       *    Importance of pheromone in ACO
       * beta : float
       *    Importance of heuristic information in ACO
       * rho : float
       *    Pheromone evaporation rate
       * w : float
       *    Inertia weight in PSO
       * c1 : float
       *    Cognitive parameter in PSO
       * c2 : float
       *    Social parameter in PSO
       * q0 : float
       *    Exploitation vs exploration balance parameter
       * max_iterations : int
       *    Maximum number of iterations
       */ this.cities = cities;
        this.n_cities = cities.length;
        this.n_ants = n_ants;
        this.n_particles = n_particles;
        this.alpha = alpha;
        this.beta = beta;
        this.rho = rho;
        this.w = w;
        this.c1 = c1;
        this.c2 = c2;
        this.q0 = q0;
        this.max_iterations = max_iterations;
        // Calculate distance matrix
        this.distance_matrix = Array(this.n_cities).fill().map(()=>Array(this.n_cities).fill(0));
        for(let i = 0; i < this.n_cities; i++){
            for(let j = 0; j < this.n_cities; j++){
                if (i !== j) {
                    this.distance_matrix[i][j] = Math.sqrt(Math.pow(cities[i][0] - cities[j][0], 2) + Math.pow(cities[i][1] - cities[j][1], 2));
                } else {
                    this.distance_matrix[i][j] = Infinity;
                }
            }
        }
        // Initialize pheromone matrix
        this.pheromone = Array(this.n_cities).fill().map(()=>Array(this.n_cities).fill(1.0));
        // Initialize heuristic information (inverse of distance)
        this.heuristic = Array(this.n_cities).fill().map(()=>Array(this.n_cities).fill(0));
        for(let i = 0; i < this.n_cities; i++){
            for(let j = 0; j < this.n_cities; j++){
                this.heuristic[i][j] = 1.0 / (this.distance_matrix[i][j] + 1e-10);
            }
        }
        // Initialize best solutions
        this.best_route = null;
        this.best_distance = Infinity;
        // Particle state
        this.particles = [];
        this.particle_velocities = [];
        this.particle_best_positions = [];
        this.particle_best_distances = Array(n_particles).fill(Infinity);
        this.global_best_position = null;
        this.global_best_distance = Infinity;
        // Initialize particles (PSO)
        this.initialize_particles();
    }
    initialize_particles() {
        /**Initialize particles with random permutations of cities*/ for(let i = 0; i < this.n_particles; i++){
            // Generate a random route
            const route = this.generateRandomPermutation(this.n_cities);
            this.particles.push(route);
            // Initialize velocity as an empty list of swaps
            this.particle_velocities.push([]);
            // Set personal best to initial position
            this.particle_best_positions.push([
                ...route
            ]);
            // Calculate initial distance
            const distance = this.calculate_route_distance(route);
            this.particle_best_distances[i] = distance;
            // Update global best if needed
            if (distance < this.global_best_distance) {
                this.global_best_distance = distance;
                this.global_best_position = [
                    ...route
                ];
            }
        }
    }
    generateRandomPermutation(n) {
        // Create array [0, 1, ..., n-1]
        const array = Array.from({
            length: n
        }, (_, i)=>i);
        // Fisher-Yates shuffle
        for(let i = array.length - 1; i > 0; i--){
            const j = Math.floor(Math.random() * (i + 1));
            [array[i], array[j]] = [
                array[j],
                array[i]
            ];
        }
        return array;
    }
    calculate_route_distance(route) {
        /**Calculate the total distance of a route*/ let distance = 0;
        for(let i = 0; i < this.n_cities; i++){
            distance += this.distance_matrix[route[i]][route[(i + 1) % this.n_cities]];
        }
        return distance;
    }
    run_aco() {
        /**Run one iteration of the ACO algorithm*/ const ant_routes = [];
        const ant_distances = [];
        // Construct solutions for each ant
        for(let ant = 0; ant < this.n_ants; ant++){
            const route = Array(this.n_cities).fill(0);
            const visited = Array(this.n_cities).fill(false);
            // Start from a random city
            const current_city = Math.floor(Math.random() * this.n_cities);
            route[0] = current_city;
            visited[current_city] = true;
            // Construct the tour
            for(let i = 1; i < this.n_cities; i++){
                let next_city;
                if (Math.random() < this.q0) {
                    // Exploitation: choose the best next city
                    const probabilities = Array(this.n_cities).fill(0);
                    for(let j = 0; j < this.n_cities; j++){
                        if (!visited[j]) {
                            probabilities[j] = Math.pow(this.pheromone[current_city][j], this.alpha) * Math.pow(this.heuristic[current_city][j], this.beta);
                        }
                    }
                    // Find city with max probability
                    let max_prob = -1;
                    let max_idx = -1;
                    for(let j = 0; j < this.n_cities; j++){
                        if (probabilities[j] > max_prob) {
                            max_prob = probabilities[j];
                            max_idx = j;
                        }
                    }
                    if (max_idx !== -1) {
                        next_city = max_idx;
                    } else {
                        // If no valid cities with pheromones, choose randomly from unvisited
                        const unvisited = [];
                        for(let j = 0; j < this.n_cities; j++){
                            if (!visited[j]) unvisited.push(j);
                        }
                        next_city = unvisited[Math.floor(Math.random() * unvisited.length)];
                    }
                } else {
                    // Exploration: choose based on probability
                    const probabilities = Array(this.n_cities).fill(0);
                    for(let j = 0; j < this.n_cities; j++){
                        if (!visited[j]) {
                            probabilities[j] = Math.pow(this.pheromone[current_city][j], this.alpha) * Math.pow(this.heuristic[current_city][j], this.beta);
                        }
                    }
                    // Normalize probabilities
                    const prob_sum = probabilities.reduce((sum, val)=>sum + val, 0);
                    if (prob_sum > 0) {
                        for(let j = 0; j < this.n_cities; j++){
                            probabilities[j] /= prob_sum;
                        }
                        // Choose next city based on probability (roulette wheel selection)
                        const r = Math.random();
                        let cum_prob = 0;
                        let selected_idx = -1;
                        for(let j = 0; j < this.n_cities; j++){
                            cum_prob += probabilities[j];
                            if (r <= cum_prob) {
                                selected_idx = j;
                                break;
                            }
                        }
                        if (selected_idx !== -1) {
                            next_city = selected_idx;
                        } else {
                            // Fallback: choose last city with non-zero probability
                            for(let j = this.n_cities - 1; j >= 0; j--){
                                if (probabilities[j] > 0) {
                                    next_city = j;
                                    break;
                                }
                            }
                        }
                    } else {
                        // If all probabilities are zero, use uniform distribution over unvisited cities
                        const unvisited = [];
                        for(let j = 0; j < this.n_cities; j++){
                            if (!visited[j]) unvisited.push(j);
                        }
                        next_city = unvisited[Math.floor(Math.random() * unvisited.length)];
                    }
                }
                route[i] = next_city;
                visited[next_city] = true;
            }
            // Calculate route distance
            const distance = this.calculate_route_distance(route);
            ant_routes.push(route);
            ant_distances.push(distance);
            // Update best solution if needed
            if (distance < this.best_distance) {
                this.best_distance = distance;
                this.best_route = [
                    ...route
                ];
            }
        }
        // Update pheromones
        // Evaporation
        for(let i = 0; i < this.n_cities; i++){
            for(let j = 0; j < this.n_cities; j++){
                this.pheromone[i][j] = (1 - this.rho) * this.pheromone[i][j];
            }
        }
        // Deposit new pheromones
        for(let ant = 0; ant < this.n_ants; ant++){
            const route = ant_routes[ant];
            const distance = ant_distances[ant];
            for(let i = 0; i < this.n_cities; i++){
                const from_city = route[i];
                const to_city = route[(i + 1) % this.n_cities];
                this.pheromone[from_city][to_city] += 1.0 / distance;
                this.pheromone[to_city][from_city] += 1.0 / distance; // Assuming symmetric TSP
            }
        }
        return [
            this.best_route,
            this.best_distance
        ];
    }
    run_pso() {
        /**Run one iteration of the PSO algorithm*/ for(let i = 0; i < this.n_particles; i++){
            // Apply velocity (swap operations)
            const route = [
                ...this.particles[i]
            ];
            // Create new velocity with inertia
            const new_velocity = [];
            if (this.particle_velocities[i].length > 0) {
                for (const swap of this.particle_velocities[i]){
                    if (Math.random() < this.w) {
                        new_velocity.push(swap);
                    }
                }
            }
            // Cognitive component (move towards personal best)
            const p_best = this.particle_best_positions[i];
            const cognitive_swaps = this.get_swaps(route, p_best);
            for (const swap of cognitive_swaps){
                if (Math.random() < this.c1 * Math.random()) {
                    new_velocity.push(swap);
                }
            }
            // Social component (move towards global best)
            const social_swaps = this.get_swaps(route, this.global_best_position);
            for (const swap of social_swaps){
                if (Math.random() < this.c2 * Math.random()) {
                    new_velocity.push(swap);
                }
            }
            // Apply velocity (swaps)
            for (const swap of new_velocity){
                const [j, k] = swap;
                [route[j], route[k]] = [
                    route[k],
                    route[j]
                ];
            }
            // Update velocity
            this.particle_velocities[i] = new_velocity;
            // Update position
            this.particles[i] = route;
            // Calculate new distance
            const distance = this.calculate_route_distance(route);
            // Update personal best
            if (distance < this.particle_best_distances[i]) {
                this.particle_best_distances[i] = distance;
                this.particle_best_positions[i] = [
                    ...route
                ];
            }
            // Update global best
            if (distance < this.global_best_distance) {
                this.global_best_distance = distance;
                this.global_best_position = [
                    ...route
                ];
            }
        }
        // Transfer information from PSO to ACO by updating pheromones based on global best
        if (this.global_best_position) {
            // Boost pheromones on global best path
            for(let i = 0; i < this.n_cities; i++){
                const from_city = this.global_best_position[i];
                const to_city = this.global_best_position[(i + 1) % this.n_cities];
                this.pheromone[from_city][to_city] += 2.0 / this.global_best_distance;
                this.pheromone[to_city][from_city] += 2.0 / this.global_best_distance;
            }
        }
        // Update best solution
        if (this.global_best_distance < this.best_distance) {
            this.best_distance = this.global_best_distance;
            this.best_route = [
                ...this.global_best_position
            ];
        }
        return [
            this.global_best_position,
            this.global_best_distance
        ];
    }
    get_swaps(source, target) {
        /**Get a list of swap operations to transform source into target*/ const swaps = [];
        const source_copy = [
            ...source
        ];
        for(let i = 0; i < this.n_cities; i++){
            if (source_copy[i] !== target[i]) {
                // Find position of target[i] in source
                const j = source_copy.indexOf(target[i]);
                // Swap
                [source_copy[i], source_copy[j]] = [
                    source_copy[j],
                    source_copy[i]
                ];
                swaps.push([
                    i,
                    j
                ]);
            }
        }
        return swaps;
    }
    optimize(verbose = true) {
        /**Run the hybrid optimization algorithm*/ const history = [];
        const iterationData = [];
        for(let iteration = 0; iteration < this.max_iterations; iteration++){
            // Run ACO
            const [aco_route, aco_distance] = this.run_aco();
            // Run PSO
            const [pso_route, pso_distance] = this.run_pso();
            // Record best distance
            history.push(this.best_distance);
            // Record detailed iteration data
            iterationData.push({
                iteration: iteration,
                aco_distance: aco_distance,
                pso_distance: pso_distance,
                best_distance: this.best_distance
            });
            if (verbose && iteration % 10 === 0) {
                console.log(`Iteration ${iteration}, Best distance: ${this.best_distance.toFixed(2)}`);
            }
        }
        // Generate route coordinates for plotting
        const routeCoordinates = this.generateRouteCoordinates(this.best_route);
        return {
            // Basic results
            route: this.best_route,
            distance: this.best_distance,
            // Data for plotting
            history: history,
            iterationData: iterationData,
            // City data
            cities: this.cities,
            cityLabels: Array.from({
                length: this.n_cities
            }, (_, i)=>i.toString()),
            // Route visualization data
            routeCoordinates: routeCoordinates,
            // Convergence plot data
            convergenceData: {
                xValues: Array.from({
                    length: this.max_iterations
                }, (_, i)=>i),
                yValues: history,
                title: 'Convergence of Hybrid ACO-PSO Algorithm',
                xLabel: 'Iteration',
                yLabel: 'Best Distance'
            },
            // Route plot data
            routePlotData: {
                cityCoordinates: this.cities,
                routePath: routeCoordinates,
                title: `Best Route - Distance: ${this.best_distance.toFixed(2)}`
            }
        };
    }
    generateRouteCoordinates(route) {
        // Generate coordinates for the route path
        const coordinates = [];
        for(let i = 0; i < route.length; i++){
            coordinates.push(this.cities[route[i]]);
        }
        // Add the first city again to close the loop
        coordinates.push(this.cities[route[0]]);
        return coordinates;
    }
}
// Example usage
function runExample() {
    // Generate random cities
    const seedRandom = (seed)=>{
        return ()=>{
            seed = (seed * 9301 + 49297) % 233280;
            return seed / 233280;
        };
    };
    const random = seedRandom(42);
    const n_cities = 20;
    const cities = Array(n_cities).fill().map(()=>[
            random() * 100,
            random() * 100
        ]);
    // Create and run the hybrid algorithm
    const hybrid = new HybridACOPSO({
        cities: cities,
        n_ants: 20,
        n_particles: 20,
        alpha: 1.0,
        beta: 5.0,
        rho: 0.5,
        w: 0.5,
        c1: 1.5,
        c2: 1.5,
        q0: 0.9,
        max_iterations: 100
    });
    const result = hybrid.optimize();
    console.log(`Best distance found: ${result.distance.toFixed(2)}`);
    console.log(`Best route: ${result.route}`);
    return result;
}
// Export the class and example function
module.exports = {
    HybridACOPSO,
    runExample
};
if (typeof globalThis.$RefreshHelpers$ === 'object' && globalThis.$RefreshHelpers !== null) {
    __turbopack_context__.k.registerExports(module, globalThis.$RefreshHelpers$);
}
}}),
"[project]/app/page.tsx [app-client] (ecmascript)": ((__turbopack_context__) => {
"use strict";

var { g: global, __dirname, k: __turbopack_refresh__, m: module } = __turbopack_context__;
{
__turbopack_context__.s({
    "default": (()=>Home)
});
var __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__ = __turbopack_context__.i("[project]/node_modules/next/dist/compiled/react/jsx-dev-runtime.js [app-client] (ecmascript)");
var __TURBOPACK__imported__module__$5b$project$5d2f$app$2f$Implementation$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__ = __turbopack_context__.i("[project]/app/Implementation.js [app-client] (ecmascript)");
var __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$index$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__ = __turbopack_context__.i("[project]/node_modules/next/dist/compiled/react/index.js [app-client] (ecmascript)");
var __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$chart$2e$js$2f$dist$2f$chart$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__$3c$locals$3e$__ = __turbopack_context__.i("[project]/node_modules/chart.js/dist/chart.js [app-client] (ecmascript) <locals>");
var __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$react$2d$chartjs$2d$2$2f$dist$2f$index$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__ = __turbopack_context__.i("[project]/node_modules/react-chartjs-2/dist/index.js [app-client] (ecmascript)");
;
var _s = __turbopack_context__.k.signature();
'use client';
;
;
;
;
__TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$chart$2e$js$2f$dist$2f$chart$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__$3c$locals$3e$__["Chart"].register(__TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$chart$2e$js$2f$dist$2f$chart$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__$3c$locals$3e$__["CategoryScale"], __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$chart$2e$js$2f$dist$2f$chart$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__$3c$locals$3e$__["LinearScale"], __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$chart$2e$js$2f$dist$2f$chart$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__$3c$locals$3e$__["PointElement"], __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$chart$2e$js$2f$dist$2f$chart$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__$3c$locals$3e$__["LineElement"], __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$chart$2e$js$2f$dist$2f$chart$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__$3c$locals$3e$__["Title"], __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$chart$2e$js$2f$dist$2f$chart$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__$3c$locals$3e$__["Tooltip"], __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$chart$2e$js$2f$dist$2f$chart$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__$3c$locals$3e$__["Legend"], __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$chart$2e$js$2f$dist$2f$chart$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__$3c$locals$3e$__["ScatterController"]);
function Home() {
    _s();
    const [cities, setCities] = (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$index$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["useState"])([]);
    const [newCityX, setNewCityX] = (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$index$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["useState"])('');
    const [newCityY, setNewCityY] = (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$index$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["useState"])('');
    const [result, setResult] = (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$index$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["useState"])(null);
    const [isOptimizing, setIsOptimizing] = (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$index$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["useState"])(false);
    const [parameters, setParameters] = (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$index$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["useState"])({
        n_ants: 20,
        n_particles: 20,
        alpha: 1.0,
        beta: 5.0,
        rho: 0.5,
        max_iterations: 100
    });
    const canvasRef = (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$index$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["useRef"])(null);
    const inputCity = ()=>{
        const x = parseFloat(newCityX);
        const y = parseFloat(newCityY);
        if (!isNaN(x) && !isNaN(y)) {
            setCities([
                ...cities,
                {
                    x,
                    y
                }
            ]);
            setNewCityX('');
            setNewCityY('');
        }
    };
    const removeCity = (index)=>{
        setCities(cities.filter((_, i)=>i !== index));
    };
    const generateRandomCities = (count)=>{
        const newCities = Array(count).fill(0).map(()=>({
                x: Math.random() * 100,
                y: Math.random() * 100
            }));
        setCities(newCities);
    };
    const clearCities = ()=>{
        setCities([]);
        setResult(null);
    };
    const runOptimization = ()=>{
        if (cities.length < 4) {
            alert("Please add at least 4 cities to run the optimization");
            return;
        }
        setIsOptimizing(true);
        // Convert cities to the format expected by HybridACOPSO
        const cityCoordinates = cities.map((city)=>[
                city.x,
                city.y
            ]);
        // Run optimization in the next tick to allow UI update
        setTimeout(()=>{
            try {
                const hybrid = new __TURBOPACK__imported__module__$5b$project$5d2f$app$2f$Implementation$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["HybridACOPSO"]({
                    cities: cityCoordinates,
                    n_ants: parameters.n_ants,
                    n_particles: parameters.n_particles,
                    alpha: parameters.alpha,
                    beta: parameters.beta,
                    rho: parameters.rho,
                    max_iterations: parameters.max_iterations
                });
                const optimizationResult = hybrid.optimize();
                setResult(optimizationResult);
            } catch (error) {
                console.error("Optimization error:", error);
                alert("An error occurred during optimization");
            } finally{
                setIsOptimizing(false);
            }
        }, 100);
    };
    // Draw cities and route on canvas
    (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$index$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["useEffect"])({
        "Home.useEffect": ()=>{
            const canvas = canvasRef.current;
            if (!canvas) return;
            const ctx = canvas.getContext('2d');
            if (!ctx) return;
            // Clear canvas
            ctx.clearRect(0, 0, canvas.width, canvas.height);
            // Scale factors to fit within canvas
            const padding = 20;
            const maxX = Math.max(...cities.map({
                "Home.useEffect.maxX": (city)=>city.x
            }["Home.useEffect.maxX"]), 100);
            const maxY = Math.max(...cities.map({
                "Home.useEffect.maxY": (city)=>city.y
            }["Home.useEffect.maxY"]), 100);
            const scaleX = (canvas.width - 2 * padding) / maxX;
            const scaleY = (canvas.height - 2 * padding) / maxY;
            // Draw cities
            cities.forEach({
                "Home.useEffect": (city, index)=>{
                    ctx.beginPath();
                    ctx.arc(padding + city.x * scaleX, padding + city.y * scaleY, 5, 0, 2 * Math.PI);
                    ctx.fillStyle = 'blue';
                    ctx.fill();
                    ctx.closePath();
                    // Add city index label
                    ctx.fillStyle = 'black';
                    ctx.font = '12px Arial';
                    ctx.fillText(`${index}`, padding + city.x * scaleX + 8, padding + city.y * scaleY + 4);
                }
            }["Home.useEffect"]);
            // Draw route if result exists and route is not null
            if (result && result.route) {
                ctx.beginPath();
                const firstCity = cities[result.route[0]];
                ctx.moveTo(padding + firstCity.x * scaleX, padding + firstCity.y * scaleY);
                // Draw route lines
                for(let i = 1; i < result.route.length; i++){
                    const city = cities[result.route[i]];
                    ctx.lineTo(padding + city.x * scaleX, padding + city.y * scaleY);
                }
                // Close the loop
                ctx.lineTo(padding + firstCity.x * scaleX, padding + firstCity.y * scaleY);
                ctx.strokeStyle = 'red';
                ctx.lineWidth = 2;
                ctx.stroke();
            }
        }
    }["Home.useEffect"], [
        cities,
        result
    ]);
    // Prepare convergence chart data
    const convergenceChartData = {
        labels: result?.convergenceData?.xValues || [],
        datasets: [
            {
                label: 'Best Distance',
                data: result?.convergenceData?.yValues || [],
                borderColor: 'rgb(75, 192, 192)',
                backgroundColor: 'rgba(75, 192, 192, 0.5)'
            }
        ]
    };
    const convergenceChartOptions = {
        responsive: true,
        plugins: {
            legend: {
                position: 'top'
            },
            title: {
                display: true,
                text: result?.convergenceData?.title || 'Convergence Chart'
            }
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
    return /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("main", {
        className: "w-full text-black min-h-screen bg-blue-100 p-8",
        children: [
            /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("h1", {
                className: "text-3xl font-serif mb-6 text-center",
                children: "Travelling Sales Man Problem using - Hybrid Ant Colony Optimization and Particle Swarm Optimization (HybridACOPSO)"
            }, void 0, false, {
                fileName: "[project]/app/page.tsx",
                lineNumber: 254,
                columnNumber: 7
            }, this),
            /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("div", {
                className: "grid grid-cols-1 md:grid-cols-2 gap-8",
                children: [
                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("div", {
                        className: "bg-gray-50 p-6 rounded-lg shadow",
                        children: [
                            /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("h2", {
                                className: "text-xl font-bold mb-4",
                                children: "Configuration"
                            }, void 0, false, {
                                fileName: "[project]/app/page.tsx",
                                lineNumber: 260,
                                columnNumber: 11
                            }, this),
                            /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("div", {
                                className: "mb-4",
                                children: [
                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("h3", {
                                        className: "font-bold mb-2",
                                        children: "Add Cities"
                                    }, void 0, false, {
                                        fileName: "[project]/app/page.tsx",
                                        lineNumber: 263,
                                        columnNumber: 13
                                    }, this),
                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("div", {
                                        className: "flex gap-2 mb-2",
                                        children: [
                                            /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("input", {
                                                type: "number",
                                                value: newCityX,
                                                onChange: (e)=>setNewCityX(e.target.value),
                                                placeholder: "X coordinate",
                                                className: "border p-2 w-32"
                                            }, void 0, false, {
                                                fileName: "[project]/app/page.tsx",
                                                lineNumber: 265,
                                                columnNumber: 15
                                            }, this),
                                            /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("input", {
                                                type: "number",
                                                value: newCityY,
                                                onChange: (e)=>setNewCityY(e.target.value),
                                                placeholder: "Y coordinate",
                                                className: "border p-2 w-32"
                                            }, void 0, false, {
                                                fileName: "[project]/app/page.tsx",
                                                lineNumber: 272,
                                                columnNumber: 15
                                            }, this),
                                            /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("button", {
                                                onClick: inputCity,
                                                className: "bg-blue-500 text-white p-2 rounded",
                                                children: "Add City"
                                            }, void 0, false, {
                                                fileName: "[project]/app/page.tsx",
                                                lineNumber: 279,
                                                columnNumber: 15
                                            }, this)
                                        ]
                                    }, void 0, true, {
                                        fileName: "[project]/app/page.tsx",
                                        lineNumber: 264,
                                        columnNumber: 13
                                    }, this),
                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("div", {
                                        className: "flex gap-2",
                                        children: [
                                            /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("button", {
                                                onClick: ()=>generateRandomCities(20),
                                                className: "bg-green-500 text-white p-2 rounded",
                                                children: "Generate 20 Random Cities"
                                            }, void 0, false, {
                                                fileName: "[project]/app/page.tsx",
                                                lineNumber: 288,
                                                columnNumber: 15
                                            }, this),
                                            /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("button", {
                                                onClick: clearCities,
                                                className: "bg-red-500 text-white p-2 rounded",
                                                children: "Clear All"
                                            }, void 0, false, {
                                                fileName: "[project]/app/page.tsx",
                                                lineNumber: 294,
                                                columnNumber: 15
                                            }, this)
                                        ]
                                    }, void 0, true, {
                                        fileName: "[project]/app/page.tsx",
                                        lineNumber: 287,
                                        columnNumber: 13
                                    }, this)
                                ]
                            }, void 0, true, {
                                fileName: "[project]/app/page.tsx",
                                lineNumber: 262,
                                columnNumber: 11
                            }, this),
                            /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("div", {
                                className: "mb-4",
                                children: [
                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("h3", {
                                        className: "font-bold mb-2",
                                        children: "Algorithm Parameters"
                                    }, void 0, false, {
                                        fileName: "[project]/app/page.tsx",
                                        lineNumber: 304,
                                        columnNumber: 13
                                    }, this),
                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("div", {
                                        className: "grid grid-cols-2 gap-2",
                                        children: [
                                            /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("div", {
                                                children: [
                                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("label", {
                                                        className: "block text-sm",
                                                        children: "Number of Ants"
                                                    }, void 0, false, {
                                                        fileName: "[project]/app/page.tsx",
                                                        lineNumber: 307,
                                                        columnNumber: 17
                                                    }, this),
                                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("input", {
                                                        type: "number",
                                                        value: parameters.n_ants,
                                                        onChange: (e)=>setParameters({
                                                                ...parameters,
                                                                n_ants: parseInt(e.target.value)
                                                            }),
                                                        className: "border p-2 w-full"
                                                    }, void 0, false, {
                                                        fileName: "[project]/app/page.tsx",
                                                        lineNumber: 308,
                                                        columnNumber: 17
                                                    }, this)
                                                ]
                                            }, void 0, true, {
                                                fileName: "[project]/app/page.tsx",
                                                lineNumber: 306,
                                                columnNumber: 15
                                            }, this),
                                            /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("div", {
                                                children: [
                                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("label", {
                                                        className: "block text-sm",
                                                        children: "Number of Particles"
                                                    }, void 0, false, {
                                                        fileName: "[project]/app/page.tsx",
                                                        lineNumber: 316,
                                                        columnNumber: 17
                                                    }, this),
                                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("input", {
                                                        type: "number",
                                                        value: parameters.n_particles,
                                                        onChange: (e)=>setParameters({
                                                                ...parameters,
                                                                n_particles: parseInt(e.target.value)
                                                            }),
                                                        className: "border p-2 w-full"
                                                    }, void 0, false, {
                                                        fileName: "[project]/app/page.tsx",
                                                        lineNumber: 317,
                                                        columnNumber: 17
                                                    }, this)
                                                ]
                                            }, void 0, true, {
                                                fileName: "[project]/app/page.tsx",
                                                lineNumber: 315,
                                                columnNumber: 15
                                            }, this),
                                            /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("div", {
                                                children: [
                                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("label", {
                                                        className: "block text-sm",
                                                        children: "Alpha (ACO)"
                                                    }, void 0, false, {
                                                        fileName: "[project]/app/page.tsx",
                                                        lineNumber: 325,
                                                        columnNumber: 17
                                                    }, this),
                                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("input", {
                                                        type: "number",
                                                        step: "0.1",
                                                        value: parameters.alpha,
                                                        onChange: (e)=>setParameters({
                                                                ...parameters,
                                                                alpha: parseFloat(e.target.value)
                                                            }),
                                                        className: "border p-2 w-full"
                                                    }, void 0, false, {
                                                        fileName: "[project]/app/page.tsx",
                                                        lineNumber: 326,
                                                        columnNumber: 17
                                                    }, this)
                                                ]
                                            }, void 0, true, {
                                                fileName: "[project]/app/page.tsx",
                                                lineNumber: 324,
                                                columnNumber: 15
                                            }, this),
                                            /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("div", {
                                                children: [
                                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("label", {
                                                        className: "block text-sm",
                                                        children: "Beta (ACO)"
                                                    }, void 0, false, {
                                                        fileName: "[project]/app/page.tsx",
                                                        lineNumber: 335,
                                                        columnNumber: 17
                                                    }, this),
                                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("input", {
                                                        type: "number",
                                                        step: "0.1",
                                                        value: parameters.beta,
                                                        onChange: (e)=>setParameters({
                                                                ...parameters,
                                                                beta: parseFloat(e.target.value)
                                                            }),
                                                        className: "border p-2 w-full"
                                                    }, void 0, false, {
                                                        fileName: "[project]/app/page.tsx",
                                                        lineNumber: 336,
                                                        columnNumber: 17
                                                    }, this)
                                                ]
                                            }, void 0, true, {
                                                fileName: "[project]/app/page.tsx",
                                                lineNumber: 334,
                                                columnNumber: 15
                                            }, this),
                                            /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("div", {
                                                children: [
                                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("label", {
                                                        className: "block text-sm",
                                                        children: "Rho (Evaporation)"
                                                    }, void 0, false, {
                                                        fileName: "[project]/app/page.tsx",
                                                        lineNumber: 345,
                                                        columnNumber: 17
                                                    }, this),
                                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("input", {
                                                        type: "number",
                                                        step: "0.1",
                                                        value: parameters.rho,
                                                        onChange: (e)=>setParameters({
                                                                ...parameters,
                                                                rho: parseFloat(e.target.value)
                                                            }),
                                                        className: "border p-2 w-full"
                                                    }, void 0, false, {
                                                        fileName: "[project]/app/page.tsx",
                                                        lineNumber: 346,
                                                        columnNumber: 17
                                                    }, this)
                                                ]
                                            }, void 0, true, {
                                                fileName: "[project]/app/page.tsx",
                                                lineNumber: 344,
                                                columnNumber: 15
                                            }, this),
                                            /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("div", {
                                                children: [
                                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("label", {
                                                        className: "block text-sm",
                                                        children: "Max Iterations"
                                                    }, void 0, false, {
                                                        fileName: "[project]/app/page.tsx",
                                                        lineNumber: 355,
                                                        columnNumber: 17
                                                    }, this),
                                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("input", {
                                                        type: "number",
                                                        value: parameters.max_iterations,
                                                        onChange: (e)=>setParameters({
                                                                ...parameters,
                                                                max_iterations: parseInt(e.target.value)
                                                            }),
                                                        className: "border p-2 w-full"
                                                    }, void 0, false, {
                                                        fileName: "[project]/app/page.tsx",
                                                        lineNumber: 356,
                                                        columnNumber: 17
                                                    }, this)
                                                ]
                                            }, void 0, true, {
                                                fileName: "[project]/app/page.tsx",
                                                lineNumber: 354,
                                                columnNumber: 15
                                            }, this)
                                        ]
                                    }, void 0, true, {
                                        fileName: "[project]/app/page.tsx",
                                        lineNumber: 305,
                                        columnNumber: 13
                                    }, this)
                                ]
                            }, void 0, true, {
                                fileName: "[project]/app/page.tsx",
                                lineNumber: 303,
                                columnNumber: 11
                            }, this),
                            /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("div", {
                                className: "mb-4",
                                children: /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("button", {
                                    onClick: runOptimization,
                                    disabled: isOptimizing || cities.length < 4,
                                    className: `w-full p-3 rounded font-bold ${isOptimizing || cities.length < 4 ? 'bg-gray-300 cursor-not-allowed' : 'bg-blue-600 text-white hover:bg-blue-700'}`,
                                    children: isOptimizing ? 'Optimizing...' : 'Run Optimization'
                                }, void 0, false, {
                                    fileName: "[project]/app/page.tsx",
                                    lineNumber: 367,
                                    columnNumber: 13
                                }, this)
                            }, void 0, false, {
                                fileName: "[project]/app/page.tsx",
                                lineNumber: 366,
                                columnNumber: 11
                            }, this),
                            /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("div", {
                                children: [
                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("h3", {
                                        className: "font-bold mb-2",
                                        children: [
                                            "City List (",
                                            cities.length,
                                            ")"
                                        ]
                                    }, void 0, true, {
                                        fileName: "[project]/app/page.tsx",
                                        lineNumber: 381,
                                        columnNumber: 13
                                    }, this),
                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("div", {
                                        className: "max-h-60 overflow-y-auto border rounded",
                                        children: /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("table", {
                                            className: "w-full",
                                            children: [
                                                /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("thead", {
                                                    className: "bg-gray-100",
                                                    children: /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("tr", {
                                                        children: [
                                                            /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("th", {
                                                                className: "p-2 text-left",
                                                                children: "Index"
                                                            }, void 0, false, {
                                                                fileName: "[project]/app/page.tsx",
                                                                lineNumber: 386,
                                                                columnNumber: 21
                                                            }, this),
                                                            /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("th", {
                                                                className: "p-2 text-left",
                                                                children: "X"
                                                            }, void 0, false, {
                                                                fileName: "[project]/app/page.tsx",
                                                                lineNumber: 387,
                                                                columnNumber: 21
                                                            }, this),
                                                            /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("th", {
                                                                className: "p-2 text-left",
                                                                children: "Y"
                                                            }, void 0, false, {
                                                                fileName: "[project]/app/page.tsx",
                                                                lineNumber: 388,
                                                                columnNumber: 21
                                                            }, this),
                                                            /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("th", {
                                                                className: "p-2 text-left",
                                                                children: "Action"
                                                            }, void 0, false, {
                                                                fileName: "[project]/app/page.tsx",
                                                                lineNumber: 389,
                                                                columnNumber: 21
                                                            }, this)
                                                        ]
                                                    }, void 0, true, {
                                                        fileName: "[project]/app/page.tsx",
                                                        lineNumber: 385,
                                                        columnNumber: 19
                                                    }, this)
                                                }, void 0, false, {
                                                    fileName: "[project]/app/page.tsx",
                                                    lineNumber: 384,
                                                    columnNumber: 17
                                                }, this),
                                                /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("tbody", {
                                                    children: [
                                                        cities.map((city, index)=>/*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("tr", {
                                                                className: "border-t",
                                                                children: [
                                                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("td", {
                                                                        className: "p-2",
                                                                        children: index
                                                                    }, void 0, false, {
                                                                        fileName: "[project]/app/page.tsx",
                                                                        lineNumber: 395,
                                                                        columnNumber: 23
                                                                    }, this),
                                                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("td", {
                                                                        className: "p-2",
                                                                        children: city.x.toFixed(2)
                                                                    }, void 0, false, {
                                                                        fileName: "[project]/app/page.tsx",
                                                                        lineNumber: 396,
                                                                        columnNumber: 23
                                                                    }, this),
                                                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("td", {
                                                                        className: "p-2",
                                                                        children: city.y.toFixed(2)
                                                                    }, void 0, false, {
                                                                        fileName: "[project]/app/page.tsx",
                                                                        lineNumber: 397,
                                                                        columnNumber: 23
                                                                    }, this),
                                                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("td", {
                                                                        className: "p-2",
                                                                        children: /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("button", {
                                                                            onClick: ()=>removeCity(index),
                                                                            className: "text-red-500 hover:text-red-700",
                                                                            children: "Remove"
                                                                        }, void 0, false, {
                                                                            fileName: "[project]/app/page.tsx",
                                                                            lineNumber: 399,
                                                                            columnNumber: 25
                                                                        }, this)
                                                                    }, void 0, false, {
                                                                        fileName: "[project]/app/page.tsx",
                                                                        lineNumber: 398,
                                                                        columnNumber: 23
                                                                    }, this)
                                                                ]
                                                            }, index, true, {
                                                                fileName: "[project]/app/page.tsx",
                                                                lineNumber: 394,
                                                                columnNumber: 21
                                                            }, this)),
                                                        cities.length === 0 && /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("tr", {
                                                            children: /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("td", {
                                                                colSpan: 4,
                                                                className: "p-2 text-center text-gray-500",
                                                                children: "No cities added yet"
                                                            }, void 0, false, {
                                                                fileName: "[project]/app/page.tsx",
                                                                lineNumber: 410,
                                                                columnNumber: 23
                                                            }, this)
                                                        }, void 0, false, {
                                                            fileName: "[project]/app/page.tsx",
                                                            lineNumber: 409,
                                                            columnNumber: 21
                                                        }, this)
                                                    ]
                                                }, void 0, true, {
                                                    fileName: "[project]/app/page.tsx",
                                                    lineNumber: 392,
                                                    columnNumber: 17
                                                }, this)
                                            ]
                                        }, void 0, true, {
                                            fileName: "[project]/app/page.tsx",
                                            lineNumber: 383,
                                            columnNumber: 15
                                        }, this)
                                    }, void 0, false, {
                                        fileName: "[project]/app/page.tsx",
                                        lineNumber: 382,
                                        columnNumber: 13
                                    }, this)
                                ]
                            }, void 0, true, {
                                fileName: "[project]/app/page.tsx",
                                lineNumber: 380,
                                columnNumber: 11
                            }, this)
                        ]
                    }, void 0, true, {
                        fileName: "[project]/app/page.tsx",
                        lineNumber: 259,
                        columnNumber: 9
                    }, this),
                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("div", {
                        children: [
                            /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("div", {
                                className: "bg-gray-50 p-6 rounded-lg shadow mb-8",
                                children: [
                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("h2", {
                                        className: "text-xl font-bold mb-4",
                                        children: "City Map & Optimal Route"
                                    }, void 0, false, {
                                        fileName: "[project]/app/page.tsx",
                                        lineNumber: 423,
                                        columnNumber: 13
                                    }, this),
                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("canvas", {
                                        ref: canvasRef,
                                        width: 500,
                                        height: 400,
                                        className: "w-full border rounded bg-white"
                                    }, void 0, false, {
                                        fileName: "[project]/app/page.tsx",
                                        lineNumber: 424,
                                        columnNumber: 13
                                    }, this),
                                    result && /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("div", {
                                        className: "mt-2 text-center font-bold text-green-600",
                                        children: [
                                            "Best Route Distance: ",
                                            result.distance.toFixed(2)
                                        ]
                                    }, void 0, true, {
                                        fileName: "[project]/app/page.tsx",
                                        lineNumber: 431,
                                        columnNumber: 15
                                    }, this)
                                ]
                            }, void 0, true, {
                                fileName: "[project]/app/page.tsx",
                                lineNumber: 422,
                                columnNumber: 11
                            }, this),
                            result && /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("div", {
                                className: "bg-gray-50 p-6 rounded-lg shadow",
                                children: [
                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("h2", {
                                        className: "text-xl font-bold mb-4",
                                        children: "Convergence Graph"
                                    }, void 0, false, {
                                        fileName: "[project]/app/page.tsx",
                                        lineNumber: 439,
                                        columnNumber: 15
                                    }, this),
                                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])(__TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$react$2d$chartjs$2d$2$2f$dist$2f$index$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["Line"], {
                                        data: convergenceChartData,
                                        options: convergenceChartOptions
                                    }, void 0, false, {
                                        fileName: "[project]/app/page.tsx",
                                        lineNumber: 440,
                                        columnNumber: 15
                                    }, this)
                                ]
                            }, void 0, true, {
                                fileName: "[project]/app/page.tsx",
                                lineNumber: 438,
                                columnNumber: 13
                            }, this)
                        ]
                    }, void 0, true, {
                        fileName: "[project]/app/page.tsx",
                        lineNumber: 421,
                        columnNumber: 9
                    }, this)
                ]
            }, void 0, true, {
                fileName: "[project]/app/page.tsx",
                lineNumber: 258,
                columnNumber: 7
            }, this),
            result && /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("div", {
                className: "mt-8 bg-gray-50 p-6 rounded-lg shadow",
                children: [
                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("h2", {
                        className: "text-xl font-bold mb-4",
                        children: "Optimization Results"
                    }, void 0, false, {
                        fileName: "[project]/app/page.tsx",
                        lineNumber: 451,
                        columnNumber: 11
                    }, this),
                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("div", {
                        className: "mb-4",
                        children: [
                            /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("h3", {
                                className: "font-bold",
                                children: "Best Route:"
                            }, void 0, false, {
                                fileName: "[project]/app/page.tsx",
                                lineNumber: 453,
                                columnNumber: 13
                            }, this),
                            /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("div", {
                                className: "p-2 bg-white border rounded overflow-x-auto",
                                children: /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("code", {
                                    children: result.route ? `${result.route.join(' → ')} → ${result.route[0]}` : 'No route found'
                                }, void 0, false, {
                                    fileName: "[project]/app/page.tsx",
                                    lineNumber: 455,
                                    columnNumber: 15
                                }, this)
                            }, void 0, false, {
                                fileName: "[project]/app/page.tsx",
                                lineNumber: 454,
                                columnNumber: 13
                            }, this)
                        ]
                    }, void 0, true, {
                        fileName: "[project]/app/page.tsx",
                        lineNumber: 452,
                        columnNumber: 11
                    }, this),
                    /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("div", {
                        children: /*#__PURE__*/ (0, __TURBOPACK__imported__module__$5b$project$5d2f$node_modules$2f$next$2f$dist$2f$compiled$2f$react$2f$jsx$2d$dev$2d$runtime$2e$js__$5b$app$2d$client$5d$__$28$ecmascript$29$__["jsxDEV"])("h3", {
                            className: "font-bold",
                            children: [
                                "Final Distance: ",
                                result.distance.toFixed(2)
                            ]
                        }, void 0, true, {
                            fileName: "[project]/app/page.tsx",
                            lineNumber: 461,
                            columnNumber: 13
                        }, this)
                    }, void 0, false, {
                        fileName: "[project]/app/page.tsx",
                        lineNumber: 460,
                        columnNumber: 11
                    }, this)
                ]
            }, void 0, true, {
                fileName: "[project]/app/page.tsx",
                lineNumber: 450,
                columnNumber: 9
            }, this)
        ]
    }, void 0, true, {
        fileName: "[project]/app/page.tsx",
        lineNumber: 253,
        columnNumber: 5
    }, this);
}
_s(Home, "q7n9/pamAlC2dIPcHB79Q1BQb50=");
_c = Home;
var _c;
__turbopack_context__.k.register(_c, "Home");
if (typeof globalThis.$RefreshHelpers$ === 'object' && globalThis.$RefreshHelpers !== null) {
    __turbopack_context__.k.registerExports(module, globalThis.$RefreshHelpers$);
}
}}),
}]);

//# sourceMappingURL=app_7392cb8d._.js.map