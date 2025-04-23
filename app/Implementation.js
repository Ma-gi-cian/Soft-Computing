class HybridACOPSO {
    constructor({
      cities,
      n_ants = 10,
      n_particles = 10,
      alpha = 1.0,
      beta = 2.0,
      rho = 0.5,
      w = 0.5,
      c1 = 1.5,
      c2 = 1.5,
      q0 = 0.9,
      max_iterations = 100
    }) {
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
       */
      this.cities = cities;
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
      this.distance_matrix = Array(this.n_cities).fill().map(() => Array(this.n_cities).fill(0));
      for (let i = 0; i < this.n_cities; i++) {
        for (let j = 0; j < this.n_cities; j++) {
          if (i !== j) {
            this.distance_matrix[i][j] = Math.sqrt(
              Math.pow(cities[i][0] - cities[j][0], 2) + 
              Math.pow(cities[i][1] - cities[j][1], 2)
            );
          } else {
            this.distance_matrix[i][j] = Infinity;
          }
        }
      }
      
      // Initialize pheromone matrix
      this.pheromone = Array(this.n_cities).fill().map(() => Array(this.n_cities).fill(1.0));
      
      // Initialize heuristic information (inverse of distance)
      this.heuristic = Array(this.n_cities).fill().map(() => Array(this.n_cities).fill(0));
      for (let i = 0; i < this.n_cities; i++) {
        for (let j = 0; j < this.n_cities; j++) {
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
      /**Initialize particles with random permutations of cities*/
      for (let i = 0; i < this.n_particles; i++) {
        // Generate a random route
        const route = this.generateRandomPermutation(this.n_cities);
        this.particles.push(route);
        
        // Initialize velocity as an empty list of swaps
        this.particle_velocities.push([]);
        
        // Set personal best to initial position
        this.particle_best_positions.push([...route]);
        
        // Calculate initial distance
        const distance = this.calculate_route_distance(route);
        this.particle_best_distances[i] = distance;
        
        // Update global best if needed
        if (distance < this.global_best_distance) {
          this.global_best_distance = distance;
          this.global_best_position = [...route];
        }
      }
    }
    
    generateRandomPermutation(n) {
      // Create array [0, 1, ..., n-1]
      const array = Array.from({ length: n }, (_, i) => i);
      
      // Fisher-Yates shuffle
      for (let i = array.length - 1; i > 0; i--) {
        const j = Math.floor(Math.random() * (i + 1));
        [array[i], array[j]] = [array[j], array[i]];
      }
      
      return array;
    }
    
    calculate_route_distance(route) {
      /**Calculate the total distance of a route*/
      let distance = 0;
      for (let i = 0; i < this.n_cities; i++) {
        distance += this.distance_matrix[route[i]][route[(i + 1) % this.n_cities]];
      }
      return distance;
    }
    
    run_aco() {
      /**Run one iteration of the ACO algorithm*/
      const ant_routes = [];
      const ant_distances = [];
      
      // Construct solutions for each ant
      for (let ant = 0; ant < this.n_ants; ant++) {
        const route = Array(this.n_cities).fill(0);
        const visited = Array(this.n_cities).fill(false);
        
        // Start from a random city
        const current_city = Math.floor(Math.random() * this.n_cities);
        route[0] = current_city;
        visited[current_city] = true;
        
        // Construct the tour
        for (let i = 1; i < this.n_cities; i++) {
          let next_city;
          
          if (Math.random() < this.q0) {
            // Exploitation: choose the best next city
            const probabilities = Array(this.n_cities).fill(0);
            for (let j = 0; j < this.n_cities; j++) {
              if (!visited[j]) {
                probabilities[j] = Math.pow(this.pheromone[current_city][j], this.alpha) * 
                                Math.pow(this.heuristic[current_city][j], this.beta);
              }
            }
            
            // Find city with max probability
            let max_prob = -1;
            let max_idx = -1;
            for (let j = 0; j < this.n_cities; j++) {
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
              for (let j = 0; j < this.n_cities; j++) {
                if (!visited[j]) unvisited.push(j);
              }
              next_city = unvisited[Math.floor(Math.random() * unvisited.length)];
            }
          } else {
            // Exploration: choose based on probability
            const probabilities = Array(this.n_cities).fill(0);
            for (let j = 0; j < this.n_cities; j++) {
              if (!visited[j]) {
                probabilities[j] = Math.pow(this.pheromone[current_city][j], this.alpha) *
                               Math.pow(this.heuristic[current_city][j], this.beta);
              }
            }
            
            // Normalize probabilities
            const prob_sum = probabilities.reduce((sum, val) => sum + val, 0);
            if (prob_sum > 0) {
              for (let j = 0; j < this.n_cities; j++) {
                probabilities[j] /= prob_sum;
              }
              
              // Choose next city based on probability (roulette wheel selection)
              const r = Math.random();
              let cum_prob = 0;
              let selected_idx = -1;
              
              for (let j = 0; j < this.n_cities; j++) {
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
                for (let j = this.n_cities - 1; j >= 0; j--) {
                  if (probabilities[j] > 0) {
                    next_city = j;
                    break;
                  }
                }
              }
            } else {
              // If all probabilities are zero, use uniform distribution over unvisited cities
              const unvisited = [];
              for (let j = 0; j < this.n_cities; j++) {
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
          this.best_route = [...route];
        }
      }
      
      // Update pheromones
      // Evaporation
      for (let i = 0; i < this.n_cities; i++) {
        for (let j = 0; j < this.n_cities; j++) {
          this.pheromone[i][j] = (1 - this.rho) * this.pheromone[i][j];
        }
      }
      
      // Deposit new pheromones
      for (let ant = 0; ant < this.n_ants; ant++) {
        const route = ant_routes[ant];
        const distance = ant_distances[ant];
        
        for (let i = 0; i < this.n_cities; i++) {
          const from_city = route[i];
          const to_city = route[(i + 1) % this.n_cities];
          this.pheromone[from_city][to_city] += 1.0 / distance;
          this.pheromone[to_city][from_city] += 1.0 / distance;  // Assuming symmetric TSP
        }
      }
      
      return [this.best_route, this.best_distance];
    }
    
    run_pso() {
      /**Run one iteration of the PSO algorithm*/
      for (let i = 0; i < this.n_particles; i++) {
        // Apply velocity (swap operations)
        const route = [...this.particles[i]];
        
        // Create new velocity with inertia
        const new_velocity = [];
        if (this.particle_velocities[i].length > 0) {
          for (const swap of this.particle_velocities[i]) {
            if (Math.random() < this.w) {
              new_velocity.push(swap);
            }
          }
        }
        
        // Cognitive component (move towards personal best)
        const p_best = this.particle_best_positions[i];
        const cognitive_swaps = this.get_swaps(route, p_best);
        for (const swap of cognitive_swaps) {
          if (Math.random() < this.c1 * Math.random()) {
            new_velocity.push(swap);
          }
        }
        
        // Social component (move towards global best)
        const social_swaps = this.get_swaps(route, this.global_best_position);
        for (const swap of social_swaps) {
          if (Math.random() < this.c2 * Math.random()) {
            new_velocity.push(swap);
          }
        }
        
        // Apply velocity (swaps)
        for (const swap of new_velocity) {
          const [j, k] = swap;
          [route[j], route[k]] = [route[k], route[j]];
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
          this.particle_best_positions[i] = [...route];
        }
        
        // Update global best
        if (distance < this.global_best_distance) {
          this.global_best_distance = distance;
          this.global_best_position = [...route];
        }
      }
      
      // Transfer information from PSO to ACO by updating pheromones based on global best
      if (this.global_best_position) {
        // Boost pheromones on global best path
        for (let i = 0; i < this.n_cities; i++) {
          const from_city = this.global_best_position[i];
          const to_city = this.global_best_position[(i + 1) % this.n_cities];
          this.pheromone[from_city][to_city] += 2.0 / this.global_best_distance;
          this.pheromone[to_city][from_city] += 2.0 / this.global_best_distance;
        }
      }
      
      // Update best solution
      if (this.global_best_distance < this.best_distance) {
        this.best_distance = this.global_best_distance;
        this.best_route = [...this.global_best_position];
      }
      
      return [this.global_best_position, this.global_best_distance];
    }
    
    get_swaps(source, target) {
      /**Get a list of swap operations to transform source into target*/
      const swaps = [];
      const source_copy = [...source];
      
      for (let i = 0; i < this.n_cities; i++) {
        if (source_copy[i] !== target[i]) {
          // Find position of target[i] in source
          const j = source_copy.indexOf(target[i]);
          // Swap
          [source_copy[i], source_copy[j]] = [source_copy[j], source_copy[i]];
          swaps.push([i, j]);
        }
      }
      
      return swaps;
    }
    
    optimize(verbose = true) {
      /**Run the hybrid optimization algorithm*/
      const history = [];
      const iterationData = [];
      
      for (let iteration = 0; iteration < this.max_iterations; iteration++) {
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
        cityLabels: Array.from({ length: this.n_cities }, (_, i) => i.toString()),
        
        // Route visualization data
        routeCoordinates: routeCoordinates,
        
        // Convergence plot data
        convergenceData: {
          xValues: Array.from({ length: this.max_iterations }, (_, i) => i),
          yValues: history,
          title: 'Convergence of Hybrid ACO-PSO Algorithm',
          xLabel: 'Iteration',
          yLabel: 'Best Distance'
        },
        
        // Route plot data
        routePlotData: {
          cityCoordinates: this.cities,
          routePath: routeCoordinates,
          title: `Best Route - Distance: ${this.best_distance.toFixed(2)}`,
        }
      };
    }
    
    generateRouteCoordinates(route) {
      // Generate coordinates for the route path
      const coordinates = [];
      for (let i = 0; i < route.length; i++) {
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
    const seedRandom = (seed) => {
      return () => {
        seed = (seed * 9301 + 49297) % 233280;
        return seed / 233280;
      };
    };
    
    const random = seedRandom(42);
    
    const n_cities = 20;
    const cities = Array(n_cities).fill().map(() => [
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
  module.exports = { HybridACOPSO, runExample };