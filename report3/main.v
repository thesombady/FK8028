import math // For math functions
import rand // For random numbers
import os // For saving to file

import vector { Vector } // Moved the vector to a separate file

// Constants needed in the simulation
const gev = 931.396 * 1e6

const speed_of_light = 299_792_458

const gram_to_ev = 5.609588845 * 1e32 / math.pow(speed_of_light * 1e10, 2)

const atom_mass = 39.948 * gev / math.pow(speed_of_light * 1e10, 2)

const density = 1.4 * gram_to_ev / 1e24

const boltzman = 8.617333262145 * 1e-5

const number_of_bins = 100


// State enum, for the state of the simulation
enum State {
	bouncing
	lattice
}

// Buffer to keep track of already computed interactions
struct Buffer {
	force Vector
	potential f64
}

// Atom struct
@[heap]
pub struct Atom {
pub mut:
	position Vector
	velocity Vector
}


// Lennard Jones potential, computes the potential and force between two atoms
fn lennard_jones(pos1 Vector, pos2 Vector, size f64) (f64, Vector) {
	r := minimum_convesion_vector(pos1, pos2, size)

	delta := math.pow(3.4 / r.norm(), 6)

	potential := 4.0 * .0104 * (delta * delta - delta)

	force := r.mul(-4.0 * 0.0104 * ( 12 * delta * delta - 6 * delta) / math.pow(r.norm(), 2))

	return potential, force
}

// Minimum convesion converts the position to the minimum value
fn minimum_convesion_vector(pos1 Vector, pos2 Vector, size f64) Vector {
	mut delta := pos2 - pos1
	delta.x += -size * math.round(delta.x / size)
	delta.y += -size * math.round(delta.y / size)
	delta.z += -size * math.round(delta.z / size)
	return delta
}

// Verlet method for the position
fn verlet_method_position(atom Atom, acceleration Vector, dt f64) Vector {
	return atom.position + atom.velocity.mul(dt) + (acceleration.mul(dt * dt/2.0))
}

// Verlet method for the velocity
fn verlet_method_velocity(atom Atom, dt f64, acc_new Vector, acc_old Vector) Vector {
	return atom.velocity + (acc_new + acc_old).mul(dt/2.0)
}

// Validate the boundary, wrapping the atoms around the box
fn validate_boundary(mut pos Vector, size f64) Vector {
	pos.x = pos.x - math.floor(pos.x / size) * size
	pos.y = pos.y - math.floor(pos.y / size) * size
	pos.z = pos.z - math.floor(pos.z / size) * size
	return pos 
}

// Initialize the atoms
// : State { Bouncing || Lattice }
// : number_of_atoms_1d { int }
fn initialize(state State, number_of_atoms_1d int) !([]Atom, f64) {
	if state == .bouncing {
		mut atoms := []Atom{}
		atoms << Atom { position: Vector.new(0.07, 0.0, 0.0), velocity: Vector.new(0.0, 0.0, 0.0) }
		atoms << Atom { position: Vector.new(6.07, 0.0, 0.0), velocity: Vector.new(0.0, 0.0, 0.0) }
		return atoms, 10.0
	}

	size := math.pow(atom_mass * math.pow(number_of_atoms_1d, 3) / density, 1.0/3.0)

	mut atoms := []Atom{}

	spacing := size / (number_of_atoms_1d - 1) * 0.8 // * 0.8

	mut average_velocity := Vector.new(0.0, 0.0, 0.0)

	for i:= 0; i < number_of_atoms_1d; i++ {
		for j:= 0; j < number_of_atoms_1d; j++ {
			for k:= 0; k < number_of_atoms_1d; k++ {
				vx := rand.f64_in_range(-1.0, 1.0)!
				vy := rand.f64_in_range(-1.0, 1.0)!
				vz := rand.f64_in_range(-1.0, 1.0)!
				pos := Vector { x: i * spacing, y: j * spacing, z: k * spacing}
				velocity := Vector { x: vx, y: vy, z: vz }.mul(1e-10) // Scaling such that the initial velocity is in Å/s
				atoms << Atom { position: pos, velocity: velocity }
				average_velocity += velocity
			}
		}
	}

	average_velocity = average_velocity.mul(1 / (math.pow(number_of_atoms_1d, 3)))

	// Drift correction
	for i in 0..atoms.len {
		atoms[i].velocity = atoms[i].velocity - average_velocity
	}

	return atoms, size
}

// Kinetic temperature
// Calculates the kinetic temperature of the system given the velocities of the atoms
fn kinetic_temperature(velocities []Vector) f64 {
	mut sum := 0.0

	for v in velocities {
		sum += v.norm() * v.norm()
	}

	return sum / (3.0 * f64(velocities.len) * boltzman) * atom_mass 
}

// Write the vector to a string
fn write_to_string(v []Vector) string {
	mut s := ''
	for j in 0..v.len {
		xyz := v[j].to_array()
		if j < v.len - 1 {
			s += '[${xyz[0]}, ${xyz[1]}, ${xyz[2]}]\t'
		} else {
			s += '[${xyz[0]}, ${xyz[1]}, ${xyz[2]}]'
		}
	}
	
	return s
}

// Save the iteration to a file (positions)
fn save_iteration(mut file os.File, v []Vector) {
	file.writeln(write_to_string(v)) or {
		panic("Could not write to file")
	}
}

// Computes the histogram for two atoms
fn compute_histogram(atom1 Atom, atom2 Atom, mut histogram []f64, dr f64, size f64) {
	r := minimum_convesion_vector(atom1.position, atom2.position, size)
	idx := int(r.norm() / dr)
	histogram[idx] += 1
}

// Divide two arrays of equal length element-wise and multipy by a factor
fn divide_array(arr1 []f64, arr2 []f64, factor f64) []f64 {
	if arr1.len != arr2.len {
		panic("The arrays must have the same length")
	}

	mut result := []f64{len: arr1.len}

	for i in 0..arr1.len {
		result[i] = factor * arr1[i] / ( arr2[i] * arr2[i] ) 
	}

	return result
}

// Save the iteration to a file (potential)
fn save_array(mut file os.File, v []f64) {
	mut s := ''
	for j in 0..v.len {
		if j < v.len - 1 {
			s += '${v[j]}\t'
		} else {
			s += '${v[j]}'
		}
	}
	file.writeln(s) or {
		panic("Could not write to file")
	}
}

// Computes the interaction between all atoms in the system.
fn compute_interactions(mut new_accelerations []Vector, mut all_potential []f64, atoms []Atom, size f64) {
	mut buffer := map[string]Buffer{}
	for i in 0..atoms.len {

		mut potential := 0.0
		mut force := Vector.new(0.0, 0.0, 0.0)

		inner: for j in 0..atoms.len {
			if i == j {
				continue inner
			}
			if '${j}:${i}' in buffer {
				potential += buffer['${j}:${i}'].potential
				force -= buffer['${j}:${i}'].force	
			} else {
				potential_, force_ := lennard_jones(atoms[i].position, atoms[j].position, size)
				potential += potential_
				force += force_
				buffer['${i}:${j}'] = Buffer {
					force: force_,
					potential: potential_
				}
			}
		}

		new_accelerations[i] += force.div(atom_mass)
		all_potential[i] = potential
	}
	buffer.clear()
}


// Run the simulation for MD with PBC
// : state { lattice || bouncing } - The state of the simulation
// : number_of_atoms_1d { int } the number of atoms in one dimension when using the lattice state.
fn run_simulation(state State, number_of_atoms_1d int)! {
	mut dt := 1e-15
	mut atoms, size := initialize(state, number_of_atoms_1d)!
	mut simulation_time := 20000

	if state == .bouncing {
		simulation_time = 5000
		dt = 5e-15
	}

	// Initialize the acceleration, position, velocity, potential and histogram for saving
	mut accelerations := []Vector{len: atoms.len}
	mut all_position := []Vector{len: atoms.len}
	mut all_velocity := []Vector{len: atoms.len}
	mut all_potential := []f64{len: atoms.len, init: 0.0}
	mut histogram := []f64{len: number_of_bins, init: 0.0}

	dr := f64(size) / f64(number_of_bins) // Casting to allow floating point division

	rs := []f64{len: number_of_bins, init: index * dr + 1e-5} // Small offset to avoid division by zero

	mut fixed_temperature := 0.0

	// Compute the initial acceleration
	compute_interactions(mut accelerations, mut all_potential, atoms, size)

	// Normalize the histogram
	for i in 0..histogram.len {
		histogram[i] /= 2.0 * f64(atoms.len)
	}

	// Compute the radial distribution function
	factor := math.pow(size, 3) / (f64(atoms.len * 4) * math.pi * dr)
	mut g := divide_array(histogram, rs, factor)

	for i in 0..atoms.len {
		all_position[i] = atoms[i].position
		all_velocity[i] = atoms[i].velocity
	}

	mut s := ''

	if state == .bouncing {
		s = 'bouncing'
	} else {
		s = 'lattice'
	}

	mut file_pos := os.create('${s}_pos.txt')!
	mut file_vel := os.create('${s}_vel.txt')!
	mut file_pot := os.create('${s}_pot.txt')!
	mut file_hist := os.create('${s}_hist.txt')!

	defer { // Defer statement to close the file towards the end
		file_pos.close()
		file_vel.close()
		file_pot.close()
		file_hist.close()
	}

	// Save the initial data, and then all the other simulated data.
	save_array(mut file_hist, g)

	// Simulation loop 
	for ts := 0; ts < simulation_time; ts ++ { 

		// Save the positions, velocities, potential and temperature to a file
		mut temp_velocity := []f64{len: atoms.len}

		for i in 0..atoms.len {
			temp_velocity[i] = math.pow(all_velocity[i].norm(), 2) * atom_mass * 0.5
		}

		save_iteration(mut file_pos, all_position)
		save_array(mut file_vel, temp_velocity)
		save_array(mut file_pot, all_potential)

		// Update position
		for i in 0..atoms.len {
			mut new_position := verlet_method_position(atoms[i], accelerations[i], dt)
			new_position = validate_boundary(mut new_position, size)
			all_position[i] = new_position 
			atoms[i].position = new_position
		}

		mut new_accelerations := []Vector{len: atoms.len}
		compute_interactions(mut new_accelerations, mut all_potential, atoms, size)

		// Update velocity
		for i in 0..atoms.len {
			all_velocity[i] = verlet_method_velocity(atoms[i], dt, new_accelerations[i], accelerations[i])	
			atoms[i].velocity = all_velocity[i] // If we rescale, we override the current velocity
		}

		// Overide the old acceleration with the new one
		accelerations = new_accelerations.clone()

		new_accelerations.clear() // clear the buffer so that we have no memory-leak

		// Rescaling the velocities to keep the temperature constant
		if ts < 2000 && ts % 50 == 0 && state == .lattice {
			temperature := kinetic_temperature(all_velocity)
			//dump(temperature)
			mut scaling_factor := math.sqrt(94.4 / temperature) // Target temperature is 94.4 K
			for j in 0..atoms.len {
				all_velocity[j] = all_velocity[j].mul(scaling_factor)
				atoms[j].velocity = all_velocity[j]
			}
		} 

		// Clear the histogram each iteration
		histogram.reset() 
	
		// Compute the histogram
		for i in 0..atoms.len {
			for j in 0..atoms.len {
				if i != j {
					compute_histogram(atoms[i], atoms[j], mut histogram, dr, size)
				}
			}
		}

		// Normalize the histogram
		for i in 0..histogram.len {
			histogram[i] /= 2.0 * f64(atoms.len)
		}

		// Compute the radial distribution function
		g = divide_array(histogram, rs, factor)

		// Save the radial distribution function to a file
		save_array(mut file_hist, g)

		println('Completed ${(f64(ts) / f64(simulation_time) * 100):.2}%')

	}
}

fn main() {
	run_simulation(State.bouncing, 5) or {
		panic("Could not run the simulation")
	}
}
