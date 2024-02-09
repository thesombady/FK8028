import math // For math functions
import rand // For random numbers
import os // For saving to file

import vector { Vector }

// Constants needed in the simulation
const gev = 931.396 * 1e6

const speed_of_light = 299_792_458

const gram_to_ev = 5.609588845 * 1e32 / math.pow(speed_of_light * 1e10, 2)

const atom_mass = 39.948 * gev / math.pow(speed_of_light * 1e10, 2)

const density = 1.4 * gram_to_ev / 1e24

const boltzman = 8.617333262145 * 1e-5


// State enum, for the state of the simulation
enum State {
	bouncing
	lattice
}

struct Placeholder {
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
	mut r := pos1 - pos2

	r.x = minimum_convesion(pos1.x, pos2.x, size)
	r.y = minimum_convesion(pos1.y, pos2.y, size)
	r.z = minimum_convesion(pos1.z, pos2.z, size)

	delta := math.pow(3.4 / r.norm(), 6)

	potential := 4.0 * .0104 * (delta * delta - delta)

	force := r.mul(-4.0 * 0.0104 * ( 12 * delta * delta - 6 * delta) / math.pow(r.norm(), 2))

	return potential, force
}

// Minimum convesion converts the position to the minimum value
fn minimum_convesion(pos1 f64, pos2 f64, size f64) f64 {
	mut delta := pos2 - pos1
	delta += -size * math.round(delta / size)
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
fn initialize(state State, number_of_atoms_1d int) ([]Atom, f64) {
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
				vx := rand.f64_in_range(-1.0, 1.0) or {
					panic("Could not generate random number")
				}
				vy := rand.f64_in_range(-1.0, 1.0) or {
					panic("Could not generate random number")
				}
				vz := rand.f64_in_range(-1.0, 1.0) or {
					panic("Could not generate random number")
				}
				pos := Vector { x: i * spacing, y: j * spacing, z: k * spacing}
				velocity := Vector { x: vx, y: vy, z: vz }.mul(1e-10)
				atoms << Atom { position: pos, velocity: velocity }
				average_velocity += velocity
			}
		}
	}

	average_velocity = average_velocity.mul(1 / (math.pow(number_of_atoms_1d, 3)))

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

fn save_iteration(mut file os.File, v []Vector) {
	file.writeln(write_to_string(v)) or {
		panic("Could not write to file")
	}
}

// Save the state of the system
fn save_state(all_position [][]Vector, all_velocity [][]Vector, all_potential [][]f64, state State) {
	mut pos_str := ''
	for i in 0..all_position.len {
		for j in 0..all_position[i].len {
			xyz := all_position[i][j].to_array()
			if j < all_position[i].len - 1 {
				pos_str += '[${xyz[0]}, ${xyz[1]}, ${xyz[2]}] \t'
			} else {
				pos_str += '[${xyz[0]}, ${xyz[1]}, ${xyz[2]}]'
			}
		}
		if i != all_position.len - 1 {
			pos_str += '\n'
		}
	}

	mut vel_str := ''
	for i in 0..all_velocity.len {
		for j in 0..all_velocity[i].len {
			xyz := all_velocity[i][j].to_array()
			if j < all_velocity[i].len - 1 {
				vel_str += '[${xyz[0]}, ${xyz[1]}, ${xyz[2]}]\t'
			} else {
				vel_str += '[${xyz[0]}, ${xyz[1]}, ${xyz[2]}]'
			}
		}
		if i != all_velocity.len - 1 {
			vel_str += ',\n'
		}
	}

	mut pot_str := '['
	for i in 0..all_potential.len {
		pot_str += '['
		for j in 0..all_potential[i].len {
			if j == all_potential[i].len - 1 {
				pot_str += '${all_potential[i][j]}]'
			} else {
				pot_str += '${all_potential[i][j]},'
			}
		}
		if i != all_potential.len - 1 {
			pot_str += ',\n'
		}
	}
	pot_str += ']'

	if state == .bouncing {

		mut file_pos := os.create('bounce_pos.txt') or {
			panic("Could not create pos file")
		}

		dump(typeof(file_pos))

		mut file_vel := os.create('bounce_vel.txt') or {
			panic("Could not create vel file")
		}

		mut file_pot := os.create('bounce_pot.txt') or {
			panic("Could not create pot file")
		}

		defer {
			file_pos.close()
			file_vel.close()
			file_pot.close()
		}

		file_pos.write(pos_str.bytes()) or {
			panic("Could not write to pos_file")
		}

		file_vel.write(vel_str.bytes()) or {
			panic("Could not write to vel_file")
		}

		file_pot.write(pot_str.bytes()) or {
			panic("Could not write to pot_file")
		}

		return
	}
	else if state == .lattice {
		mut file_pos := os.create('lattice_pos.txt') or {
			panic("Could not create pos file")
		}

		mut file_vel := os.create('lattice_vel.txt') or {
			panic("Could not create vel file")
		}

		mut file_pot := os.create('lattice_pot.txt') or {
			panic("Could not create pot file")
		}

		defer {
			file_pos.close()
			file_vel.close()
			file_pot.close()
		}

		file_pos.write(pos_str.bytes()) or {
			panic("Could not write to pos_file")
		}

		file_vel.write(vel_str.bytes()) or {
			panic("Could not write to vel_file")
		}

		file_pot.write(pot_str.bytes()) or {
			panic("Could not write to pot_file")
		}	
	}
}

// Run the simulation for MD, with random velocities and PBC.
fn run_simulation(state State, number_of_atoms_1d int) {

	mut dt := 1e-15

	mut atoms, size := initialize(state, number_of_atoms_1d)

	mut simulation_time := 5000

	if state == .bouncing {
		simulation_time = 10000
		dt = 5e-15
	}

	mut accelerations := []Vector{len: atoms.len}

	mut all_position := [][]Vector{len: simulation_time, init: []Vector{len: atoms.len}}
	mut all_velocity := [][]Vector{len: simulation_time, init: []Vector{len: atoms.len}}
	mut all_potential := [][]f64{len: simulation_time, init: []f64{len: atoms.len, init: 0.0}}
	
	// Compute the initial acceleration
	for i in 0..atoms.len {
		mut potential := 0.0
		mut force := Vector.new(0.0, 0.0, 0.0)
		for j in 0..atoms.len {
			if i != j {
				potential_, force_ := lennard_jones(atoms[i].position, atoms[j].position, size)
				potential += potential_
				force += force_
			}
		}

		accelerations[i] += force.div(atom_mass)

		all_potential[0][i] = potential
		
	}

	// Simulation loop 
	for ts := 0; ts < simulation_time; ts ++ { 

		t := ts * dt

		// Update position
		for i in 0..atoms.len {
   
			mut new_position := verlet_method_position(atoms[i], accelerations[i], dt)

			new_position = validate_boundary(mut new_position, size)

			all_position[ts][i] = new_position 

			atoms[i].position = new_position
		}

		// Compute new acceleration

		mut computed := map[string]Placeholder{}

		mut new_accelerations := []Vector{len: atoms.len}

		for i in 0..atoms.len {

			mut potential := 0.0
			mut force := Vector.new(0.0, 0.0, 0.0)

			inner: for j in 0..atoms.len {
				if i == j {
					continue inner
				}
				if '${j}:${i}' in computed {
					potential -= computed['${j}:${i}'].potential
					force -= computed['${j}:${i}'].force	
				} else {
					potential_, force_ := lennard_jones(atoms[i].position, atoms[j].position, size)
					potential += potential_
					force += force_
					computed['${i}:${j}'] = Placeholder {
						force: force_,
						potential: potential_
					}
				}
			}

			new_accelerations[i] += force.div(atom_mass)
			all_potential[ts][i] = potential
		}

		computed.clear()

		// Update velocity
		
		for i in 0..atoms.len {
			all_velocity[ts][i] = verlet_method_velocity(atoms[i], dt, new_accelerations[i], accelerations[i])	
		}

		// Rescaling the velocities to keep the temperature constant
		if ts < 2000 && ts % 50 == 0 {
			temperature := kinetic_temperature(all_velocity[ts])
			//dump(temperature)
			mut scaling_factor := math.sqrt(94.4 / temperature)
			for j in 0..atoms.len {
				all_velocity[ts][j] = all_velocity[ts][j].mul(scaling_factor)
				atoms[j].velocity = all_velocity[ts][j]
			}
		} else {
			dump(kinetic_temperature(all_velocity[ts]))
			for j in 0..atoms.len {
				atoms[j].velocity = all_velocity[ts][j]
			}
		}


		// Overide the old acceleration with the new one

		accelerations = new_accelerations.clone()

		new_accelerations.clear() // clear the buffer so that we have no memory-leak

		println('Completed ${(f64(ts) / f64(simulation_time) * 100):.2}%')

	}

	// Saving the positions, velocities, and potential to a file format for visualization

	save_state(all_position, all_velocity, all_potential, state)
}

fn main() {
	run_simulation(State.bouncing, 5)
}
