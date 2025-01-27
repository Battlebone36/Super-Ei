from code.algorithms.genetic import Genetic
from code.classes.protein import *
import copy

class Gen_1(Genetic):
    def run(self, store_step_stability: bool=False, population_size: int = 30, mutation_probability: float = 0.001, nominator: int = 50, verbose: bool = False) -> Protein:
        """
        Genetic algorithm that mimics natural selection to find the optimal folded protein.
        """
        if verbose:
            print("Starting algorithm")

        # Initialize variables
        population: dict[int, tuple[Protein, list[int]]] = {}
        generations: int = 1000
        iteration_limit = self.max_iterations

        # Create initial population - generation 0
        population, protein_id = self.create_initial_population(population_size)
        if verbose:
            print("initial population created")

        # Track best solution
        best_protein = min(population.values(), key=lambda p: p[0].stability())
        generations_without_change = 0
        store_max_protein = self.protein


        for generation in range(generations):
            if self.iterations >= iteration_limit:
                break
            
            if verbose:
                print(f"Generation {generation + 1}/{generations}")

            # Create new population
            new_population: dict[int, tuple[Protein, list[int]]] = {}
            
            # Fill the new population with offsprings
            while len(new_population) < population_size:
                # Select the "best" two proteins from the population based on stability
                parent1 = self.tournament_selection(population, nominator)
                parent2 = self.tournament_selection(population, nominator)
                
                self.iterations += 2
                if self.iterations >= iteration_limit:
                    break
                if store_step_stability:
                    self.store_steps_stability()
                    
                # Create the offsprings with crossover
                child1_folds, child2_folds = self.build_offsprings(parent1, parent2)
            
                # Apply mutations
                child1_folds = self.mutation_on_sequence(child1_folds, mutation_probability)
                child2_folds = self.mutation_on_sequence(child2_folds, mutation_probability)

                # Apply the folds to the offsprings
                child1 = Protein(self.protein.sequence)
                child2 = Protein(self.protein.sequence)

                folded_child1 = self.fold_protein_by_sequence(child1, child1_folds)
                folded_child2 = self.fold_protein_by_sequence(child2, child2_folds)

                # Add offsprings to the new population
                new_population[protein_id] = (folded_child1, child1_folds)
                protein_id += 1

                new_population[protein_id] = (folded_child2, child2_folds)
                protein_id += 1
                temporarily_best_protein = min(new_population.values(), key=lambda p: p[0].stability())

                # Track solution for data store
                current_best_stability = store_max_protein.stability()
                new_stability = temporarily_best_protein[0].stability()
                if current_best_stability > new_stability:
                    self.protein = copy.deepcopy(temporarily_best_protein[0])
                
            population = new_population

            # Update best solution
            current_best_protein = min(population.values(), key=lambda p: p[0].stability())
            if current_best_protein[0].stability() < best_protein[0].stability():
                best_protein = current_best_protein
                generations_without_change = 0
            else:
                generations_without_change += 1

            if verbose:
                print(f"Best protein stability at generation {generation + 1}: {best_protein[0].stability()}")

            # Early exit if there is no improvement for consecutive generations

            if generations_without_change > 100 and not store_step_stability:
                if verbose:
                    print(f"No improvement for 100 generations. Early stopping at generation {generation + 1}.")
                break
        
        if verbose:
            print(f"Final best protein stability: {best_protein[0].stability()}")

        self.protein = best_protein[0]
        return best_protein[0]
    
class Gen_2(Genetic):
    def run(self, store_step_stability: bool=False, population_size: int = 30, mutation_probability: float = 0.01, nominator: int = 70, verbose: bool = False) -> Protein:
        """
        Genetic algorithm that mimics natural selection to find the optimal folded protein.
        """
        if verbose:
            print("Starting algorithm")

        # Initialize variables
        population: dict[int, tuple[Protein, list[int]]] = {}
        generations: int = 1000
        iteration_limit = self.max_iterations

        # Create initial population - generation 0
        population, protein_id = self.create_initial_population(population_size)
        if verbose:
            print("initial population created")

        # Track best solution
        best_protein = min(population.values(), key=lambda p: p[0].stability())
        generations_without_change = 0
        store_max_protein = self.protein


        for generation in range(generations):
            if self.iterations >= iteration_limit:
                break
            
            if verbose:
                print(f"Generation {generation + 1}/{generations}")

            # Create new population
            new_population: dict[int, tuple[Protein, list[int]]] = {}
            
            # Fill the new population with offsprings
            while len(new_population) < population_size:
                # Select the "best" two proteins from the population based on stability
                parent1 = self.tournament_selection(population, nominator)
                parent2 = self.tournament_selection(population, nominator)
                
                self.iterations += 2
                if self.iterations >= iteration_limit:
                    break
                if store_step_stability:
                    self.store_steps_stability()
                    
                # Create the offsprings with crossover
                child1_folds, child2_folds = self.build_offsprings(parent1, parent2)
            
                # Apply mutations
                child1_folds = self.mutation_on_sequence(child1_folds, mutation_probability)
                child2_folds = self.mutation_on_sequence(child2_folds, mutation_probability)

                # Apply the folds to the offsprings
                child1 = Protein(self.protein.sequence)
                child2 = Protein(self.protein.sequence)

                folded_child1 = self.fold_protein_by_sequence(child1, child1_folds)
                folded_child2 = self.fold_protein_by_sequence(child2, child2_folds)

                # Add offsprings to the new population
                new_population[protein_id] = (folded_child1, child1_folds)
                protein_id += 1

                new_population[protein_id] = (folded_child2, child2_folds)
                protein_id += 1
                temporarily_best_protein = min(new_population.values(), key=lambda p: p[0].stability())

                # Track solution for data store
                current_best_stability = store_max_protein.stability()
                new_stability = temporarily_best_protein[0].stability()
                if current_best_stability > new_stability:
                    self.protein = copy.deepcopy(temporarily_best_protein[0])
                
            population = new_population

            # Update best solution
            current_best_protein = min(population.values(), key=lambda p: p[0].stability())
            if current_best_protein[0].stability() < best_protein[0].stability():
                best_protein = current_best_protein
                generations_without_change = 0
            else:
                generations_without_change += 1

            if verbose:
                print(f"Best protein stability at generation {generation + 1}: {best_protein[0].stability()}")

            # Early exit if there is no improvement for consecutive generations

            if generations_without_change > 100 and not store_step_stability:
                if verbose:
                    print(f"No improvement for 100 generations. Early stopping at generation {generation + 1}.")
                break
        
        if verbose:
            print(f"Final best protein stability: {best_protein[0].stability()}")

        self.protein = best_protein[0]
        return best_protein[0]
    
class Gen_3(Genetic):
    def run(self, store_step_stability: bool=False, population_size: int = 50, mutation_probability: float = 0.01, nominator: int = 40, verbose: bool = False) -> Protein:
        """
        Genetic algorithm that mimics natural selection to find the optimal folded protein.
        """
        if verbose:
            print("Starting algorithm")

        # Initialize variables
        population: dict[int, tuple[Protein, list[int]]] = {}
        generations: int = 1000
        iteration_limit = self.max_iterations

        # Create initial population - generation 0
        population, protein_id = self.create_initial_population(population_size)
        if verbose:
            print("initial population created")

        # Track best solution
        best_protein = min(population.values(), key=lambda p: p[0].stability())
        generations_without_change = 0
        store_max_protein = self.protein


        for generation in range(generations):
            if self.iterations >= iteration_limit:
                break
            
            if verbose:
                print(f"Generation {generation + 1}/{generations}")

            # Create new population
            new_population: dict[int, tuple[Protein, list[int]]] = {}
            
            # Fill the new population with offsprings
            while len(new_population) < population_size:
                # Select the "best" two proteins from the population based on stability
                parent1 = self.tournament_selection(population, nominator)
                parent2 = self.tournament_selection(population, nominator)
                
                self.iterations += 2
                if self.iterations >= iteration_limit:
                    break
                if store_step_stability:
                    self.store_steps_stability()
                    
                # Create the offsprings with crossover
                child1_folds, child2_folds = self.build_offsprings(parent1, parent2)
            
                # Apply mutations
                child1_folds = self.mutation_on_sequence(child1_folds, mutation_probability)
                child2_folds = self.mutation_on_sequence(child2_folds, mutation_probability)

                # Apply the folds to the offsprings
                child1 = Protein(self.protein.sequence)
                child2 = Protein(self.protein.sequence)

                folded_child1 = self.fold_protein_by_sequence(child1, child1_folds)
                folded_child2 = self.fold_protein_by_sequence(child2, child2_folds)

                # Add offsprings to the new population
                new_population[protein_id] = (folded_child1, child1_folds)
                protein_id += 1

                new_population[protein_id] = (folded_child2, child2_folds)
                protein_id += 1
                temporarily_best_protein = min(new_population.values(), key=lambda p: p[0].stability())

                # Track solution for data store
                current_best_stability = store_max_protein.stability()
                new_stability = temporarily_best_protein[0].stability()
                if current_best_stability > new_stability:
                    self.protein = copy.deepcopy(temporarily_best_protein[0])
                
            population = new_population

            # Update best solution
            current_best_protein = min(population.values(), key=lambda p: p[0].stability())
            if current_best_protein[0].stability() < best_protein[0].stability():
                best_protein = current_best_protein
                generations_without_change = 0
            else:
                generations_without_change += 1

            if verbose:
                print(f"Best protein stability at generation {generation + 1}: {best_protein[0].stability()}")

            # Early exit if there is no improvement for consecutive generations

            if generations_without_change > 100 and not store_step_stability:
                if verbose:
                    print(f"No improvement for 100 generations. Early stopping at generation {generation + 1}.")
                break
        
        if verbose:
            print(f"Final best protein stability: {best_protein[0].stability()}")

        self.protein = best_protein[0]
        return best_protein[0]
    
class Gen_4(Genetic):
    def run(self, store_step_stability: bool=False, population_size: int = 50, mutation_probability: float = 0.01, nominator: int = 50, verbose: bool = False) -> Protein:
        """
        Genetic algorithm that mimics natural selection to find the optimal folded protein.
        """
        if verbose:
            print("Starting algorithm")

        # Initialize variables
        population: dict[int, tuple[Protein, list[int]]] = {}
        generations: int = 1000
        iteration_limit = self.max_iterations

        # Create initial population - generation 0
        population, protein_id = self.create_initial_population(population_size)
        if verbose:
            print("initial population created")

        # Track best solution
        best_protein = min(population.values(), key=lambda p: p[0].stability())
        generations_without_change = 0
        store_max_protein = self.protein


        for generation in range(generations):
            if self.iterations >= iteration_limit:
                break
            
            if verbose:
                print(f"Generation {generation + 1}/{generations}")

            # Create new population
            new_population: dict[int, tuple[Protein, list[int]]] = {}
            
            # Fill the new population with offsprings
            while len(new_population) < population_size:
                # Select the "best" two proteins from the population based on stability
                parent1 = self.tournament_selection(population, nominator)
                parent2 = self.tournament_selection(population, nominator)
                
                self.iterations += 2
                if self.iterations >= iteration_limit:
                    break
                if store_step_stability:
                    self.store_steps_stability()
                    
                # Create the offsprings with crossover
                child1_folds, child2_folds = self.build_offsprings(parent1, parent2)
            
                # Apply mutations
                child1_folds = self.mutation_on_sequence(child1_folds, mutation_probability)
                child2_folds = self.mutation_on_sequence(child2_folds, mutation_probability)

                # Apply the folds to the offsprings
                child1 = Protein(self.protein.sequence)
                child2 = Protein(self.protein.sequence)

                folded_child1 = self.fold_protein_by_sequence(child1, child1_folds)
                folded_child2 = self.fold_protein_by_sequence(child2, child2_folds)

                # Add offsprings to the new population
                new_population[protein_id] = (folded_child1, child1_folds)
                protein_id += 1

                new_population[protein_id] = (folded_child2, child2_folds)
                protein_id += 1
                temporarily_best_protein = min(new_population.values(), key=lambda p: p[0].stability())

                # Track solution for data store
                current_best_stability = store_max_protein.stability()
                new_stability = temporarily_best_protein[0].stability()
                if current_best_stability > new_stability:
                    self.protein = copy.deepcopy(temporarily_best_protein[0])
                
            population = new_population

            # Update best solution
            current_best_protein = min(population.values(), key=lambda p: p[0].stability())
            if current_best_protein[0].stability() < best_protein[0].stability():
                best_protein = current_best_protein
                generations_without_change = 0
            else:
                generations_without_change += 1

            if verbose:
                print(f"Best protein stability at generation {generation + 1}: {best_protein[0].stability()}")

            # Early exit if there is no improvement for consecutive generations

            if generations_without_change > 100 and not store_step_stability:
                if verbose:
                    print(f"No improvement for 100 generations. Early stopping at generation {generation + 1}.")
                break
        
        if verbose:
            print(f"Final best protein stability: {best_protein[0].stability()}")

        self.protein = best_protein[0]
        return best_protein[0]
    
class Gen_5(Genetic):
    def run(self, store_step_stability: bool=False, population_size: int = 70, mutation_probability: float = 0.01, nominator: int = 60, verbose: bool = False) -> Protein:
        """
        Genetic algorithm that mimics natural selection to find the optimal folded protein.
        """
        if verbose:
            print("Starting algorithm")

        # Initialize variables
        population: dict[int, tuple[Protein, list[int]]] = {}
        generations: int = 1000
        iteration_limit = self.max_iterations

        # Create initial population - generation 0
        population, protein_id = self.create_initial_population(population_size)
        if verbose:
            print("initial population created")

        # Track best solution
        best_protein = min(population.values(), key=lambda p: p[0].stability())
        generations_without_change = 0
        store_max_protein = self.protein


        for generation in range(generations):
            if self.iterations >= iteration_limit:
                break
            
            if verbose:
                print(f"Generation {generation + 1}/{generations}")

            # Create new population
            new_population: dict[int, tuple[Protein, list[int]]] = {}
            
            # Fill the new population with offsprings
            while len(new_population) < population_size:
                # Select the "best" two proteins from the population based on stability
                parent1 = self.tournament_selection(population, nominator)
                parent2 = self.tournament_selection(population, nominator)
                
                self.iterations += 2
                if self.iterations >= iteration_limit:
                    break
                if store_step_stability:
                    self.store_steps_stability()
                    
                # Create the offsprings with crossover
                child1_folds, child2_folds = self.build_offsprings(parent1, parent2)
            
                # Apply mutations
                child1_folds = self.mutation_on_sequence(child1_folds, mutation_probability)
                child2_folds = self.mutation_on_sequence(child2_folds, mutation_probability)

                # Apply the folds to the offsprings
                child1 = Protein(self.protein.sequence)
                child2 = Protein(self.protein.sequence)

                folded_child1 = self.fold_protein_by_sequence(child1, child1_folds)
                folded_child2 = self.fold_protein_by_sequence(child2, child2_folds)

                # Add offsprings to the new population
                new_population[protein_id] = (folded_child1, child1_folds)
                protein_id += 1

                new_population[protein_id] = (folded_child2, child2_folds)
                protein_id += 1
                temporarily_best_protein = min(new_population.values(), key=lambda p: p[0].stability())

                # Track solution for data store
                current_best_stability = store_max_protein.stability()
                new_stability = temporarily_best_protein[0].stability()
                if current_best_stability > new_stability:
                    self.protein = copy.deepcopy(temporarily_best_protein[0])
                
            population = new_population

            # Update best solution
            current_best_protein = min(population.values(), key=lambda p: p[0].stability())
            if current_best_protein[0].stability() < best_protein[0].stability():
                best_protein = current_best_protein
                generations_without_change = 0
            else:
                generations_without_change += 1

            if verbose:
                print(f"Best protein stability at generation {generation + 1}: {best_protein[0].stability()}")

            # Early exit if there is no improvement for consecutive generations

            if generations_without_change > 100 and not store_step_stability:
                if verbose:
                    print(f"No improvement for 100 generations. Early stopping at generation {generation + 1}.")
                break
        
        if verbose:
            print(f"Final best protein stability: {best_protein[0].stability()}")

        self.protein = best_protein[0]
        return best_protein[0]
    
class Gen_6(Genetic):
    def run(self, store_step_stability: bool=False, population_size: int = 70, mutation_probability: float = 0.01, nominator: int = 30, verbose: bool = False) -> Protein:
        """
        Genetic algorithm that mimics natural selection to find the optimal folded protein.
        """
        if verbose:
            print("Starting algorithm")

        # Initialize variables
        population: dict[int, tuple[Protein, list[int]]] = {}
        generations: int = 1000
        iteration_limit = self.max_iterations

        # Create initial population - generation 0
        population, protein_id = self.create_initial_population(population_size)
        if verbose:
            print("initial population created")

        # Track best solution
        best_protein = min(population.values(), key=lambda p: p[0].stability())
        generations_without_change = 0
        store_max_protein = self.protein


        for generation in range(generations):
            if self.iterations >= iteration_limit:
                break
            
            if verbose:
                print(f"Generation {generation + 1}/{generations}")

            # Create new population
            new_population: dict[int, tuple[Protein, list[int]]] = {}
            
            # Fill the new population with offsprings
            while len(new_population) < population_size:
                # Select the "best" two proteins from the population based on stability
                parent1 = self.tournament_selection(population, nominator)
                parent2 = self.tournament_selection(population, nominator)
                
                self.iterations += 2
                if self.iterations >= iteration_limit:
                    break
                if store_step_stability:
                    self.store_steps_stability()
                    
                # Create the offsprings with crossover
                child1_folds, child2_folds = self.build_offsprings(parent1, parent2)
            
                # Apply mutations
                child1_folds = self.mutation_on_sequence(child1_folds, mutation_probability)
                child2_folds = self.mutation_on_sequence(child2_folds, mutation_probability)

                # Apply the folds to the offsprings
                child1 = Protein(self.protein.sequence)
                child2 = Protein(self.protein.sequence)

                folded_child1 = self.fold_protein_by_sequence(child1, child1_folds)
                folded_child2 = self.fold_protein_by_sequence(child2, child2_folds)

                # Add offsprings to the new population
                new_population[protein_id] = (folded_child1, child1_folds)
                protein_id += 1

                new_population[protein_id] = (folded_child2, child2_folds)
                protein_id += 1
                temporarily_best_protein = min(new_population.values(), key=lambda p: p[0].stability())

                # Track solution for data store
                current_best_stability = store_max_protein.stability()
                new_stability = temporarily_best_protein[0].stability()
                if current_best_stability > new_stability:
                    self.protein = copy.deepcopy(temporarily_best_protein[0])
                
            population = new_population

            # Update best solution
            current_best_protein = min(population.values(), key=lambda p: p[0].stability())
            if current_best_protein[0].stability() < best_protein[0].stability():
                best_protein = current_best_protein
                generations_without_change = 0
            else:
                generations_without_change += 1

            if verbose:
                print(f"Best protein stability at generation {generation + 1}: {best_protein[0].stability()}")

            # Early exit if there is no improvement for consecutive generations

            if generations_without_change > 100 and not store_step_stability:
                if verbose:
                    print(f"No improvement for 100 generations. Early stopping at generation {generation + 1}.")
                break
        
        if verbose:
            print(f"Final best protein stability: {best_protein[0].stability()}")

        self.protein = best_protein[0]
        return best_protein[0]
    
class Gen_7(Genetic):
    def run(self, store_step_stability: bool=False, population_size: int = 70, mutation_probability: float = 0.01, nominator: int = 70, verbose: bool = False) -> Protein:
        """
        Genetic algorithm that mimics natural selection to find the optimal folded protein.
        """
        if verbose:
            print("Starting algorithm")

        # Initialize variables
        population: dict[int, tuple[Protein, list[int]]] = {}
        generations: int = 1000
        iteration_limit = self.max_iterations

        # Create initial population - generation 0
        population, protein_id = self.create_initial_population(population_size)
        if verbose:
            print("initial population created")

        # Track best solution
        best_protein = min(population.values(), key=lambda p: p[0].stability())
        generations_without_change = 0
        store_max_protein = self.protein


        for generation in range(generations):
            if self.iterations >= iteration_limit:
                break
            
            if verbose:
                print(f"Generation {generation + 1}/{generations}")

            # Create new population
            new_population: dict[int, tuple[Protein, list[int]]] = {}
            
            # Fill the new population with offsprings
            while len(new_population) < population_size:
                # Select the "best" two proteins from the population based on stability
                parent1 = self.tournament_selection(population, nominator)
                parent2 = self.tournament_selection(population, nominator)
                
                self.iterations += 2
                if self.iterations >= iteration_limit:
                    break
                if store_step_stability:
                    self.store_steps_stability()
                    
                # Create the offsprings with crossover
                child1_folds, child2_folds = self.build_offsprings(parent1, parent2)
            
                # Apply mutations
                child1_folds = self.mutation_on_sequence(child1_folds, mutation_probability)
                child2_folds = self.mutation_on_sequence(child2_folds, mutation_probability)

                # Apply the folds to the offsprings
                child1 = Protein(self.protein.sequence)
                child2 = Protein(self.protein.sequence)

                folded_child1 = self.fold_protein_by_sequence(child1, child1_folds)
                folded_child2 = self.fold_protein_by_sequence(child2, child2_folds)

                # Add offsprings to the new population
                new_population[protein_id] = (folded_child1, child1_folds)
                protein_id += 1

                new_population[protein_id] = (folded_child2, child2_folds)
                protein_id += 1
                temporarily_best_protein = min(new_population.values(), key=lambda p: p[0].stability())

                # Track solution for data store
                current_best_stability = store_max_protein.stability()
                new_stability = temporarily_best_protein[0].stability()
                if current_best_stability > new_stability:
                    self.protein = copy.deepcopy(temporarily_best_protein[0])
                
            population = new_population

            # Update best solution
            current_best_protein = min(population.values(), key=lambda p: p[0].stability())
            if current_best_protein[0].stability() < best_protein[0].stability():
                best_protein = current_best_protein
                generations_without_change = 0
            else:
                generations_without_change += 1

            if verbose:
                print(f"Best protein stability at generation {generation + 1}: {best_protein[0].stability()}")

            # Early exit if there is no improvement for consecutive generations

            if generations_without_change > 100 and not store_step_stability:
                if verbose:
                    print(f"No improvement for 100 generations. Early stopping at generation {generation + 1}.")
                break
        
        if verbose:
            print(f"Final best protein stability: {best_protein[0].stability()}")

        self.protein = best_protein[0]
        return best_protein[0]
    
class Gen_8(Genetic):
    def run(self, store_step_stability: bool=False, population_size: int = 90, mutation_probability: float = 0.01, nominator: int = 30, verbose: bool = False) -> Protein:
        """
        Genetic algorithm that mimics natural selection to find the optimal folded protein.
        """
        if verbose:
            print("Starting algorithm")

        # Initialize variables
        population: dict[int, tuple[Protein, list[int]]] = {}
        generations: int = 1000
        iteration_limit = self.max_iterations

        # Create initial population - generation 0
        population, protein_id = self.create_initial_population(population_size)
        if verbose:
            print("initial population created")

        # Track best solution
        best_protein = min(population.values(), key=lambda p: p[0].stability())
        generations_without_change = 0
        store_max_protein = self.protein


        for generation in range(generations):
            if self.iterations >= iteration_limit:
                break
            
            if verbose:
                print(f"Generation {generation + 1}/{generations}")

            # Create new population
            new_population: dict[int, tuple[Protein, list[int]]] = {}
            
            # Fill the new population with offsprings
            while len(new_population) < population_size:
                # Select the "best" two proteins from the population based on stability
                parent1 = self.tournament_selection(population, nominator)
                parent2 = self.tournament_selection(population, nominator)
                
                self.iterations += 2
                if self.iterations >= iteration_limit:
                    break
                if store_step_stability:
                    self.store_steps_stability()
                    
                # Create the offsprings with crossover
                child1_folds, child2_folds = self.build_offsprings(parent1, parent2)
            
                # Apply mutations
                child1_folds = self.mutation_on_sequence(child1_folds, mutation_probability)
                child2_folds = self.mutation_on_sequence(child2_folds, mutation_probability)

                # Apply the folds to the offsprings
                child1 = Protein(self.protein.sequence)
                child2 = Protein(self.protein.sequence)

                folded_child1 = self.fold_protein_by_sequence(child1, child1_folds)
                folded_child2 = self.fold_protein_by_sequence(child2, child2_folds)

                # Add offsprings to the new population
                new_population[protein_id] = (folded_child1, child1_folds)
                protein_id += 1

                new_population[protein_id] = (folded_child2, child2_folds)
                protein_id += 1
                temporarily_best_protein = min(new_population.values(), key=lambda p: p[0].stability())

                # Track solution for data store
                current_best_stability = store_max_protein.stability()
                new_stability = temporarily_best_protein[0].stability()
                if current_best_stability > new_stability:
                    self.protein = copy.deepcopy(temporarily_best_protein[0])
                
            population = new_population

            # Update best solution
            current_best_protein = min(population.values(), key=lambda p: p[0].stability())
            if current_best_protein[0].stability() < best_protein[0].stability():
                best_protein = current_best_protein
                generations_without_change = 0
            else:
                generations_without_change += 1

            if verbose:
                print(f"Best protein stability at generation {generation + 1}: {best_protein[0].stability()}")

            # Early exit if there is no improvement for consecutive generations

            if generations_without_change > 100 and not store_step_stability:
                if verbose:
                    print(f"No improvement for 100 generations. Early stopping at generation {generation + 1}.")
                break
        
        if verbose:
            print(f"Final best protein stability: {best_protein[0].stability()}")

        self.protein = best_protein[0]
        return best_protein[0]
    
class Gen_9(Genetic):
    def run(self, store_step_stability: bool=False, population_size: int = 90, mutation_probability: float = 0.01, nominator: int = 70, verbose: bool = False) -> Protein:
        """
        Genetic algorithm that mimics natural selection to find the optimal folded protein.
        """
        if verbose:
            print("Starting algorithm")

        # Initialize variables
        population: dict[int, tuple[Protein, list[int]]] = {}
        generations: int = 1000
        iteration_limit = self.max_iterations

        # Create initial population - generation 0
        population, protein_id = self.create_initial_population(population_size)
        if verbose:
            print("initial population created")

        # Track best solution
        best_protein = min(population.values(), key=lambda p: p[0].stability())
        generations_without_change = 0
        store_max_protein = self.protein


        for generation in range(generations):
            if self.iterations >= iteration_limit:
                break
            
            if verbose:
                print(f"Generation {generation + 1}/{generations}")

            # Create new population
            new_population: dict[int, tuple[Protein, list[int]]] = {}
            
            # Fill the new population with offsprings
            while len(new_population) < population_size:
                # Select the "best" two proteins from the population based on stability
                parent1 = self.tournament_selection(population, nominator)
                parent2 = self.tournament_selection(population, nominator)
                
                self.iterations += 2
                if self.iterations >= iteration_limit:
                    break
                if store_step_stability:
                    self.store_steps_stability()
                    
                # Create the offsprings with crossover
                child1_folds, child2_folds = self.build_offsprings(parent1, parent2)
            
                # Apply mutations
                child1_folds = self.mutation_on_sequence(child1_folds, mutation_probability)
                child2_folds = self.mutation_on_sequence(child2_folds, mutation_probability)

                # Apply the folds to the offsprings
                child1 = Protein(self.protein.sequence)
                child2 = Protein(self.protein.sequence)

                folded_child1 = self.fold_protein_by_sequence(child1, child1_folds)
                folded_child2 = self.fold_protein_by_sequence(child2, child2_folds)

                # Add offsprings to the new population
                new_population[protein_id] = (folded_child1, child1_folds)
                protein_id += 1

                new_population[protein_id] = (folded_child2, child2_folds)
                protein_id += 1
                temporarily_best_protein = min(new_population.values(), key=lambda p: p[0].stability())

                # Track solution for data store
                current_best_stability = store_max_protein.stability()
                new_stability = temporarily_best_protein[0].stability()
                if current_best_stability > new_stability:
                    self.protein = copy.deepcopy(temporarily_best_protein[0])
                
            population = new_population

            # Update best solution
            current_best_protein = min(population.values(), key=lambda p: p[0].stability())
            if current_best_protein[0].stability() < best_protein[0].stability():
                best_protein = current_best_protein
                generations_without_change = 0
            else:
                generations_without_change += 1

            if verbose:
                print(f"Best protein stability at generation {generation + 1}: {best_protein[0].stability()}")

            # Early exit if there is no improvement for consecutive generations

            if generations_without_change > 100 and not store_step_stability:
                if verbose:
                    print(f"No improvement for 100 generations. Early stopping at generation {generation + 1}.")
                break
        
        if verbose:
            print(f"Final best protein stability: {best_protein[0].stability()}")

        self.protein = best_protein[0]
        return best_protein[0]
 
        """
        Genetic algorithm that mimics natural selection to find the optimal folded protein.
        """
        if verbose:
            print("Starting algorithm")

        # Initialize variables
        population: dict[int, tuple[Protein, list[int]]] = {}
        generations: int = 1000
        iteration_limit = self.max_iterations

        # Create initial population - generation 0
        population, protein_id = self.create_initial_population(population_size)
        if verbose:
            print("initial population created")

        # Track best solution
        best_protein = min(population.values(), key=lambda p: p[0].stability())
        generations_without_change = 0
        store_max_protein = self.protein


        for generation in range(generations):
            if self.iterations >= iteration_limit:
                break
            
            if verbose:
                print(f"Generation {generation + 1}/{generations}")

            # Create new population
            new_population: dict[int, tuple[Protein, list[int]]] = {}
            
            # Fill the new population with offsprings
            while len(new_population) < population_size:
                # Select the "best" two proteins from the population based on stability
                parent1 = self.tournament_selection(population, nominator)
                parent2 = self.tournament_selection(population, nominator)
                
                self.iterations += 2
                if self.iterations >= iteration_limit:
                    break
                if store_step_stability:
                    self.store_steps_stability()
                    
                # Create the offsprings with crossover
                child1_folds, child2_folds = self.build_offsprings(parent1, parent2)
            
                # Apply mutations
                child1_folds = self.mutation_on_sequence(child1_folds, mutation_probability)
                child2_folds = self.mutation_on_sequence(child2_folds, mutation_probability)

                # Apply the folds to the offsprings
                child1 = Protein(self.protein.sequence)
                child2 = Protein(self.protein.sequence)

                folded_child1 = self.fold_protein_by_sequence(child1, child1_folds)
                folded_child2 = self.fold_protein_by_sequence(child2, child2_folds)

                # Add offsprings to the new population
                new_population[protein_id] = (folded_child1, child1_folds)
                protein_id += 1

                new_population[protein_id] = (folded_child2, child2_folds)
                protein_id += 1
                temporarily_best_protein = min(new_population.values(), key=lambda p: p[0].stability())

                # Track solution for data store
                current_best_stability = store_max_protein.stability()
                new_stability = temporarily_best_protein[0].stability()
                if current_best_stability > new_stability:
                    self.protein = copy.deepcopy(temporarily_best_protein[0])
                
            population = new_population

            # Update best solution
            current_best_protein = min(population.values(), key=lambda p: p[0].stability())
            if current_best_protein[0].stability() < best_protein[0].stability():
                best_protein = current_best_protein
                generations_without_change = 0
            else:
                generations_without_change += 1

            if verbose:
                print(f"Best protein stability at generation {generation + 1}: {best_protein[0].stability()}")

            # Early exit if there is no improvement for consecutive generations

            if generations_without_change > 100 and not store_step_stability:
                if verbose:
                    print(f"No improvement for 100 generations. Early stopping at generation {generation + 1}.")
                break
        
        if verbose:
            print(f"Final best protein stability: {best_protein[0].stability()}")

        self.protein = best_protein[0]
        return best_protein[0]