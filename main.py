import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches

class HybridACOPSO:
    def __init__(self, cities, n_ants=10, n_particles=10, alpha=1.0, beta=2.0, rho=0.5, 
                 w=0.5, c1=1.5, c2=1.5, q0=0.9, max_iterations=100):
        """
        Initialize the hybrid ACO-PSO algorithm.
        
        Parameters:
        -----------
        cities : numpy.ndarray
            Coordinates of cities as (x, y) pairs
        n_ants : int
            Number of ants in the colony
        n_particles : int
            Number of particles in the swarm
        alpha : float
            Importance of pheromone in ACO
        beta : float
            Importance of heuristic information in ACO
        rho : float
            Pheromone evaporation rate
        w : float
            Inertia weight in PSO
        c1 : float
            Cognitive parameter in PSO
        c2 : float
            Social parameter in PSO
        q0 : float
            Exploitation vs exploration balance parameter
        max_iterations : int
            Maximum number of iterations
        """
        self.cities = cities
        self.n_cities = len(cities)
        self.n_ants = n_ants
        self.n_particles = n_particles
        self.alpha = alpha
        self.beta = beta
        self.rho = rho
        self.w = w
        self.c1 = c1
        self.c2 = c2
        self.q0 = q0
        self.max_iterations = max_iterations
        
        # Calculate distance matrix
        self.distance_matrix = np.zeros((self.n_cities, self.n_cities))
        for i in range(self.n_cities):
            for j in range(self.n_cities):
                if i != j:
                    self.distance_matrix[i, j] = np.sqrt(
                        np.sum((cities[i] - cities[j]) ** 2)
                    )
                else:
                    self.distance_matrix[i, j] = np.inf
        
        # Initialize pheromone matrix
        self.pheromone = np.ones((self.n_cities, self.n_cities))
        
        # Initialize heuristic information (inverse of distance)
        self.heuristic = 1.0 / (self.distance_matrix + 1e-10)
        
        # Initialize best solutions
        self.best_route = None
        self.best_distance = float('inf')
        
        # Particle state
        self.particles = []
        self.particle_velocities = []
        self.particle_best_positions = []
        self.particle_best_distances = np.full(n_particles, float('inf'))
        self.global_best_position = None
        self.global_best_distance = float('inf')
        
        # Initialize particles (PSO)
        self.initialize_particles()
    
    def initialize_particles(self):
        """Initialize particles with random permutations of cities"""
        for i in range(self.n_particles):
            # Generate a random route
            route = np.random.permutation(self.n_cities)
            self.particles.append(route)
            
            # Initialize velocity as an empty list of swaps
            self.particle_velocities.append([])
            
            # Set personal best to initial position
            self.particle_best_positions.append(route.copy())
            
            # Calculate initial distance
            distance = self.calculate_route_distance(route)
            self.particle_best_distances[i] = distance
            
            # Update global best if needed
            if distance < self.global_best_distance:
                self.global_best_distance = distance
                self.global_best_position = route.copy()
    
    def calculate_route_distance(self, route):
        """Calculate the total distance of a route"""
        distance = 0
        for i in range(self.n_cities):
            distance += self.distance_matrix[route[i], route[(i + 1) % self.n_cities]]
        return distance
    
    def run_aco(self):
        """Run one iteration of the ACO algorithm"""
        ant_routes = []
        ant_distances = []
        
        # Construct solutions for each ant
        for ant in range(self.n_ants):
            route = np.zeros(self.n_cities, dtype=int)
            visited = np.zeros(self.n_cities, dtype=bool)
            
            # Start from a random city
            current_city = np.random.randint(0, self.n_cities)
            route[0] = current_city
            visited[current_city] = True
            
            # Construct the tour
            for i in range(1, self.n_cities):
                if np.random.random() < self.q0:
                    # Exploitation: choose the best next city
                    probabilities = np.zeros(self.n_cities)
                    for j in range(self.n_cities):
                        if not visited[j]:
                            probabilities[j] = (self.pheromone[current_city, j] ** self.alpha) * \
                                             (self.heuristic[current_city, j] ** self.beta)
                    
                    # Prevent division by zero by checking if any valid cities
                    if np.sum(probabilities) > 0:
                        next_city = np.argmax(probabilities)
                    else:
                        # If no valid cities with pheromones, choose randomly from unvisited
                        unvisited = np.where(~visited)[0]
                        next_city = np.random.choice(unvisited)
                else:
                    # Exploration: choose based on probability
                    probabilities = np.zeros(self.n_cities)
                    for j in range(self.n_cities):
                        if not visited[j]:
                            probabilities[j] = (self.pheromone[current_city, j] ** self.alpha) * \
                                             (self.heuristic[current_city, j] ** self.beta)
                    
                    # Normalize probabilities
                    prob_sum = np.sum(probabilities)
                    if prob_sum > 0:
                        probabilities = probabilities / prob_sum
                    else:
                        # If all probabilities are zero, use uniform distribution over unvisited cities
                        unvisited = np.where(~visited)[0]
                        for j in unvisited:
                            probabilities[j] = 1.0 / len(unvisited)
                    
                    # Make sure probabilities sum to 1 (avoiding floating point issues)
                    probabilities = probabilities / np.sum(probabilities)
                    
                    # Choose next city based on probability
                    valid_indices = np.where(probabilities > 0)[0]
                    probabilities_valid = probabilities[valid_indices]
                    
                    # Final safety check to ensure valid probability distribution
                    if len(valid_indices) > 0 and np.isclose(np.sum(probabilities_valid), 1.0):
                        next_city = np.random.choice(valid_indices, p=probabilities_valid)
                    else:
                        # Fallback: choose randomly from unvisited
                        unvisited = np.where(~visited)[0]
                        next_city = np.random.choice(unvisited)
                
                route[i] = next_city
                visited[next_city] = True
                current_city = next_city
            
            # Calculate route distance
            distance = self.calculate_route_distance(route)
            ant_routes.append(route)
            ant_distances.append(distance)
            
            # Update best solution if needed
            if distance < self.best_distance:
                self.best_distance = distance
                self.best_route = route.copy()
        
        # Update pheromones
        # Evaporation
        self.pheromone = (1 - self.rho) * self.pheromone
        
        # Deposit new pheromones
        for ant in range(self.n_ants):
            route = ant_routes[ant]
            distance = ant_distances[ant]
            
            for i in range(self.n_cities):
                from_city = route[i]
                to_city = route[(i + 1) % self.n_cities]
                self.pheromone[from_city, to_city] += 1.0 / distance
                self.pheromone[to_city, from_city] += 1.0 / distance  # Assuming symmetric TSP
        
        return self.best_route, self.best_distance
    
    def run_pso(self):
        """Run one iteration of the PSO algorithm"""
        for i in range(self.n_particles):
            # Apply velocity (swap operations)
            route = self.particles[i].copy()
            
            # Create new velocity with inertia
            new_velocity = []
            if len(self.particle_velocities[i]) > 0:
                for swap in self.particle_velocities[i]:
                    if np.random.random() < self.w:
                        new_velocity.append(swap)
            
            # Cognitive component (move towards personal best)
            p_best = self.particle_best_positions[i]
            cognitive_swaps = self.get_swaps(route, p_best)
            for swap in cognitive_swaps:
                if np.random.random() < self.c1 * np.random.random():
                    new_velocity.append(swap)
            
            # Social component (move towards global best)
            social_swaps = self.get_swaps(route, self.global_best_position)
            for swap in social_swaps:
                if np.random.random() < self.c2 * np.random.random():
                    new_velocity.append(swap)
            
            # Apply velocity (swaps)
            for swap in new_velocity:
                j, k = swap
                route[j], route[k] = route[k], route[j]
            
            # Update velocity
            self.particle_velocities[i] = new_velocity
            
            # Update position
            self.particles[i] = route
            
            # Calculate new distance
            distance = self.calculate_route_distance(route)
            
            # Update personal best
            if distance < self.particle_best_distances[i]:
                self.particle_best_distances[i] = distance
                self.particle_best_positions[i] = route.copy()
            
            # Update global best
            if distance < self.global_best_distance:
                self.global_best_distance = distance
                self.global_best_position = route.copy()
        
        # Transfer information from PSO to ACO by updating pheromones based on global best
        if self.global_best_position is not None:
            # Boost pheromones on global best path
            for i in range(self.n_cities):
                from_city = self.global_best_position[i]
                to_city = self.global_best_position[(i + 1) % self.n_cities]
                self.pheromone[from_city, to_city] += 2.0 / self.global_best_distance
                self.pheromone[to_city, from_city] += 2.0 / self.global_best_distance
        
        # Update best solution
        if self.global_best_distance < self.best_distance:
            self.best_distance = self.global_best_distance
            self.best_route = self.global_best_position.copy()
        
        return self.global_best_position, self.global_best_distance
    
    def get_swaps(self, source, target):
        """Get a list of swap operations to transform source into target"""
        swaps = []
        source_copy = source.copy()
        
        for i in range(self.n_cities):
            if source_copy[i] != target[i]:
                # Find position of target[i] in source
                j = np.where(source_copy == target[i])[0][0]
                # Swap
                source_copy[i], source_copy[j] = source_copy[j], source_copy[i]
                swaps.append((i, j))
        
        return swaps
    
    def optimize(self, verbose=True):
        """Run the hybrid optimization algorithm"""
        history = []
        
        for iteration in range(self.max_iterations):
            # Run ACO
            aco_route, aco_distance = self.run_aco()
            
            # Run PSO
            pso_route, pso_distance = self.run_pso()
            
            # Record best distance
            history.append(self.best_distance)
            
            if verbose and iteration % 10 == 0:
                print(f"Iteration {iteration}, Best distance: {self.best_distance:.2f}")
        
        return self.best_route, self.best_distance, history
    
    def plot_route(self, route=None, title="Best Route"):
        """Plot the route on a 2D plane"""
        if route is None:
            route = self.best_route
        
        plt.figure(figsize=(10, 6))
        
        # Plot cities
        plt.scatter(self.cities[:, 0], self.cities[:, 1], s=100, c='blue')
        
        # Plot route
        route_coords = self.cities[route]
        route_coords = np.append(route_coords, [route_coords[0]], axis=0)  # Close the loop
        
        plt.plot(route_coords[:, 0], route_coords[:, 1], 'r-', alpha=0.7)
        
        # Annotate cities
        for i, (x, y) in enumerate(self.cities):
            plt.annotate(str(i), (x, y), xytext=(10, 10), textcoords='offset points')
        
        plt.title(f"{title} - Distance: {self.calculate_route_distance(route):.2f}")
        plt.grid(True)
        plt.show()
    
    def plot_convergence(self, history):
        """Plot the convergence of the algorithm"""
        plt.figure(figsize=(10, 6))
        plt.plot(history)
        plt.xlabel('Iteration')
        plt.ylabel('Best Distance')
        plt.title('Convergence of Hybrid ACO-PSO Algorithm')
        plt.grid(True)
        plt.show()


# Example usage
if __name__ == "__main__":
    # Generate random cities
    np.random.seed(42)
    n_cities = 20
    cities = np.random.rand(n_cities, 2) * 100
    
    # Create and run the hybrid algorithm
    hybrid = HybridACOPSO(
        cities=cities,
        n_ants=20,
        n_particles=20,
        alpha=1.0,
        beta=5.0,
        rho=0.5,
        w=0.5,
        c1=1.5,
        c2=1.5,
        q0=0.9,
        max_iterations=100
    )
    
    best_route, best_distance, history = hybrid.optimize()
    
    print(f"Best distance found: {best_distance:.2f}")
    print(f"Best route: {best_route}")
    
    # Plot results
    hybrid.plot_route()
    hybrid.plot_convergence(history)
