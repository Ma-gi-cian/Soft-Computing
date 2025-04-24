# Hybrid ACO-PSO Algorithm for the Traveling Salesman Problem



## Overview

This repo contains the ``main.py``  file containing the code towards the approach of combining both **Ant Colony Optimization** and  **Swarm Particle Optimization** . The aim towards this approach is to achieve faster and more efficient convergence via using the combined approach of both the algorithms.

1. __Ant colony Optimization__ : Allows for local optimization
2. __Swarm Particle Optimization__ : Allows for global optimization
 

## Algorithm Description

### Hybrid ACO-PSO

The implementation combines two powerful meta heuristic approach :

1. **Ant Colony Optimization (ACO)**: A population-based approach inspired by the foraging behavior of ants, which use pheromone trails to communicate and find optimal paths.

2. **Particle Swarm Optimization (PSO)**: An optimization technique inspired by social behavior of bird flocking or fish schooling, where particles move in the search space guided by their own experience and the swarm's collective knowledge.

The hybrid algorithm runs both optimization techniques in alternating fashion, with information sharing between them to enhance exploration and exploitation capabilities.

## Features

- Solves the Traveling Salesman Problem (TSP)
- Combines ACO and PSO approaches for improved performance
- Configurable parameters for both algorithms
- Supports cities with 2D coordinates

## Implementation Details

### `HybridACOPSO` Class

#### Initialization Parameters

- `cities`: 2D numpy array containing city coordinates
- `n_ants`: Number of ants in the colony (default: 10)
- `n_particles`: Number of particles in the swarm (default: 10)
- `alpha`: Importance of pheromone in ACO (default: 1.0)
- `beta`: Importance of heuristic information (inverse of distance) in ACO (default: 2.0)
- `rho`: Pheromone evaporation rate (default: 0.5)
- `w`: Inertia weight in PSO (default: 0.5)
- `c1`: Cognitive parameter in PSO (default: 1.5)
- `c2`: Social parameter in PSO (default: 1.5)
- `q0`: Exploitation vs exploration balance parameter (default: 0.9)
- `max_iterations`: Maximum number of iterations (default: 100)

#### Key Components

1. **Initialization**:
   - Creates distance matrix between all cities
   - Initializes pheromone matrix with ones
   - Computes heuristic information as inverse of distances
   - Initializes PSO particles with random permutations of cities

2. **ACO Implementation**:
   - Ants construct solutions by probabilistically selecting the next city based on pheromone levels and heuristic information
   - Uses a q0 parameter to balance between exploitation (choosing the best next city) and exploration (probabilistic choice)
   - Updates pheromone levels after all ants complete their tours
   - Includes pheromone evaporation to prevent premature convergence

3. **PSO Implementation**:
   - Represents solutions as permutations of cities
   - Velocities are stored as sequences of swap operations
   - Updates particle positions by applying a series of swaps
   - Includes cognitive component (personal best) and social component (global best)
   - Transfers information to ACO by boosting pheromone levels on the global best path

4. **Optimization Process**:
   - Alternates between ACO and PSO for the specified number of iterations
   - Tracks the best solution found and convergence history
   - Provides visualization of the best route and convergence curve

### Key Methods

- `initialize_particles()`: Creates initial random permutations for particles
- `calculate_route_distance()`: Computes the total distance of a route
- `run_aco()`: Executes one iteration of the ACO algorithm
- `run_pso()`: Executes one iteration of the PSO algorithm
- `get_swaps()`: Identifies swap operations to transform one permutation into another
- `optimize()`: Runs the full optimization process, alternating between ACO and PSO
- `plot_route()`: Visualizes the best route found
- `plot_convergence()`: Displays the convergence history

## Algorithm Workflow

1. **Initialization**:
   - Calculate distances between all cities
   - Initialize pheromone levels and PSO particles

2. **Main Loop**:
   - For each iteration:
     - Run one iteration of ACO
     - Run one iteration of PSO
     - Update best solution if improved
     - Record convergence history

3. **ACO Procedure**:
   - For each ant:
     - Start at a random city
     - Construct a tour by selecting cities based on pheromone and distance
     - Update best solution if improved
   - Update pheromone levels:
     - Evaporate existing pheromones
     - Deposit new pheromones based on ants' tour quality

4. **PSO Procedure**:
   - For each particle:
     - Apply velocity (swap operations)
     - Update position
     - Update personal best and global best if improved
   - Transfer information to ACO by boosting pheromones on global best path

5. **Information Exchange**:
   - PSO influences ACO through pheromone updates
   - ACO influences PSO through the exchange of good solutions

## Example Usage

```python
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
```

## Parameter Tuning

For optimal performance, consider adjusting these parameters:

- **ACO parameters**:
  - `alpha`: Higher values give more importance to pheromone trails (exploitation)
  - `beta`: Higher values give more importance to shorter distances (greedy behavior)
  - `rho`: Controls pheromone evaporation rate
  - `q0`: Higher values promote exploitation over exploration

- **PSO parameters**:
  - `w`: Inertia weight (controls influence of previous velocity)
  - `c1`: Cognitive parameter (attraction to personal best)
  - `c2`: Social parameter (attraction to global best)

## Visualization

The implementation includes two visualization functions:

1. `plot_route()`: Displays the best route found, with cities as blue dots and the route as red lines
2. `plot_convergence()`: Shows the convergence curve (best distance vs. iteration)


## Code 
```python
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
```

1.  The fitness function is demonstrated via the function named `calculate_route_distance(route)` . This function computes the total distance of the tour, which both the algorithms PSO and ACO try to minimize. 
```python
def calculate_route_distance(self, route):
	distance = 0
	for i in range(self.n_cities):
		distance += self.distance_matrix[route[i], route[(i+1)%self.n_cities]]
	return distance 
## Fitness = total length of the route ( goal is to minimize it )
```

2. __Declaration Phase__
- `Distance Matrix: Initialize a distance matrix with all zeros. Compute all the distances. If I is not J then calcualte the euclidean distance or else put np.inf ( infinity out there )  .`
- `Initialize the Pheromone matrix with np.ones (all ones ) `
- `self.heuristic = inverse of the distance matrix ( from ant colony optimization ) `
- `Initialize the best routes using self.best_route = None and self.best_distance = infinity `

3. ### __Initialize Particles__
- `Initialize the particles  `
- `route = np.random.permutation(self.n_cities) -> generates a random route. np.random.permutation just shuffles all the values inside of the given list. `
- `Initialize an empty particle_velocities, particle_best_positions.`
- `Calculate the distance of the route i.e., the permutatied route that is generated using the np.random.permutation`
- `Update the personal best distance and the global best distance for the particles.`

4.  **Route Distance Calculation** 
-  `Takes route as an input. and calculates the distane of the entire route.`

5. ### __Run ACO__
- `Initialize the ant_routes [integer type] and the ant_distances [ boolean type ] `
- `start from a random city  set the first city inside route matrix to current city and set visited at that city to true`
- `The code inside the for loop it iterates over all the present cities except the current selected city. `
- `if np.random.random() - random value is greater than self.q0 ( q0 is used to establish the balance between exploration and exploitation ) : this line is used to make a probabilistic decision. `
- `An array called probabilities is initialized with zeros. This is used to store the attractiveness towards the next city to be choosen. `
- `Now the J loop that iterated through all the possible cities.`
- `check if it is not visited then calculate the probabilities of  j = ( self.pheromone[current_city, j] ** self.alpha * (self.heuristic[current_city, j ]**beta `
- `Here alpha is the importance on the pheromone and beta is the importance on the heuristic `
- `check if the probability is not equal to zero ( could or could not happen dont know but during debugging it was said to make this change ) `
- `If by any chance there is no attractive city to be choosen then we choose randomly via a unifrom distribution function to select a city`
- `The next city is choosen using  np.random.choice(valid_indices, p=probabilities_valid). valid_indices = the cities that can be visited and p = the probabilites that are assosiated with those cities, and np.random.choice will give the city (its an inbuilt function in numpy ).`
- `current_city = next_city`
- `calcualte the distance between the cities and add it.`
- `self.pheromone = ( 1 - self.rho ) * self.pheromone // factor by which the pheromone evalporates.`
- `We deposite the pheromone. We increase the better pheromone via add 1/distance `

6. ### __PSO__

- `Iterate through all the particles `
- `initialize an empty velocity array.`
- `for swap in self.particle_velocities[i] we iterate through all the particles velocities and check if np.random.random() < self.w `
- `w is the inertial weight that controles the influence of the previous weight of the particle`
- `we calculate the p_best that is the personal best of the particles and then the global_best of the entire system. This method is similar to how the pso algorithm works`

> <span style="background-color: #e6f7ff; padding: 4px; border-radius: 4px;">V<sub>new</sub> = V<sub>old</sub> + C1 * rand() * ( pBest - CurrentValue ) + C2 * rand() * ( gBest - CurrentValue )</span>

> `C1 = Cognitive Parameter { Influence of the personal best on the velocity calculation } `
`C2 = Social Parameter { Influence of the global best on the velocity calculation }`

- `The above formula is used for it `
- `This interacts via that if the pso algorithm has found a global best then it will boost the pheromone between the cities. If checks for the global best and if it has a global best from its own side then it will boost the pheromone between those cities. `
- `Similar to the ACO self.pheromone[from_city, to_city] += 2.0 / self.global_best_distance `
- `The ACO did it via 1.0 / probabilities and the PSO does it via the self.global_best_distance `
- `We make these two algorithms interact via make the trails more attractive to the ants using the pso's best found distance. `

