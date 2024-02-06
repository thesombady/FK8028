import math // For math functions
import rand // For random numbers
import os // For saving to file

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

// Vector struct, for 3D vectors
pub struct Vector {
pub mut:
	x f64
	y f64
	z f64
}

pub fn Vector.new(x f64, y f64, z f64) Vector {
	return Vector { x: x, y: y, z: z }
}

@[inline]
pub fn (v Vector) + (u Vector) Vector {
	return Vector { x: v.x + u.x, y: v.y + u.y, z: v.z + u.z }
}


@[inline]
pub fn (v Vector) - (u Vector) Vector {
	return Vector { x: v.x - u.x, y: v.y - u.y, z: v.z - u.z }
}

pub fn (v Vector) mul(s f64) Vector {
	return Vector { x: v.x * s, y: v.y * s, z: v.z * s }
}

pub fn (v Vector) norm() f64 {
	return math.sqrt(v.x * v.x + v.y * v.y + v.z * v.z)
}

pub fn (v Vector) div(s f64) Vector {
	if s == 0.0 {
		panic("Division by zero")
	}

	return Vector { x: v.x / s, y: v.y / s, z: v.z / s }
}

pub fn (v Vector) get(i int) ?f64 {
	return match i {
		0 { v.x }
		1 { v.y }
		2 { v.z }
		else { none }
	}
}

pub fn (v Vector) to_array() []f64 {
	return [v.x, v.y, v.z]
}

// Atom struct
@[heap]
pub struct Atom {
pub:
	name string
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
	return  atom.position + atom.velocity + (acceleration.div(2)).mul(math.pow(dt, 2))
}

// Verlet method for the velocity
fn verlet_method_velocity(atom Atom, dt f64, acc_new Vector, acc_old Vector) Vector {
	return atom.velocity + ((acc_new + acc_old).div(2)).mul(dt)
}

// Validate the boundary, wrapping the atoms around the box
fn validate_boundary(mut pos Vector, size f64) Vector {

	pos.x = pos.x - math.floor(pos.x / size) * size
	pos.y = pos.y - math.floor(pos.y / size) * size
	pos.z = pos.z - math.floor(pos.z / size) * size

	return pos 
}

// Initialize the atoms
fn initialize(state State, number_of_atoms_1d int) ([]Atom, f64) {
	if state == .bouncing {
		mut atoms := []Atom{}
		atoms << Atom { name: "Ar:0", position: Vector.new(0.07, 0.0, 0.0), velocity: Vector.new(0.0, 0.0, 0.0) }
		atoms << Atom { name: "Ar:1", position: Vector.new(6.07, 0.0, 0.0), velocity: Vector.new(0.0, 0.0, 0.0) }
		return atoms, 10.0
	}


	size := math.pow(atom_mass * math.pow(number_of_atoms_1d, 3) / density, 1.0/3.0)

	dump(size)

	mut atoms := []Atom{}

	spacing := size / (number_of_atoms_1d - 1) * 0.8 // * 0.8

	mut counter := 0
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
				atoms << Atom { name: "Ar:${counter}", position: pos, velocity: velocity }
				counter += 1
			}
		}
	}
	return atoms, size
}

// Save to csv: Save the position, velocity and potential energy of the atoms to a csv file
fn save_to_csv(all_atoms map[string]map[string]Vector, state State, z int) {
	if state == .bouncing {
		mut file := os.create('iteration/bt${z}.csv') or {
			panic("Could not create file")
		}

		mut file_string := 'Atom, X, Y, Z, Vx, Vy, Vz, Pot\n'
		for name, _ in all_atoms {
			position := (all_atoms[name]["position"]).to_array()
			velocity := (all_atoms[name]["velocity"]).to_array()
			potential := ((all_atoms[name]["potential"]).to_array())[0]
			file_string += '${name}, ${position[0]}, ${position[1]}, ${position[2]}, ${velocity[0]}, ${velocity[1]}, ${velocity[2]}, ${potential}\n'
		}

		file.write(file_string.bytes()) or {
			panic("Could not write to file")
		}
		return 
	}	

	mut file := os.create('iteration/lattice${z}.csv') or {
		panic("Could not create file")
	}
	
	mut file_string := 'Atom, X, Y, Z, Vx, Vy, Vz, Pot\n'

	for name, _ in all_atoms {
		position := (all_atoms[name]["position"]).to_array()
		velocity := (all_atoms[name]["velocity"]).to_array()
		potential := ((all_atoms[name]["potential"]).to_array())[0]
		file_string += '${name}, ${position[0]}, ${position[1]}, ${position[2]}, ${velocity[0]}, ${velocity[1]}, ${velocity[2]}, ${potential}\n'
	}

	file.write(file_string.bytes()) or {
		panic("Could not write to file")
	}
}

// Kinetic temperature
// Calculates the kinetic temperature of the system given the velocities of the atoms
fn kinetic_temperature(velocities []Vector) f64 {
	mut sum := 0.0
	
	for v in velocities {
		sum += v.norm() * v.norm()
	}

	return sum / (3 * velocities.len) * atom_mass * boltzman
}

// Run the simulation for MD, with random velocities and PBC.
fn run_simulation(state State, number_of_atoms_1d int) {
	mut simulaton_time := 1e-12

	mut dt := 1e-15

	mut atoms, size := initialize(state, number_of_atoms_1d)

	if state == .bouncing {
		simulaton_time = 1e-11	
		dt = 0.8e-15
	}

	mut t := 0.0

	mut accelerations := []Vector{len: atoms.len}

	mut all_atoms := map[string]map[string]Vector{}

	for atom in atoms {
		all_atoms[atom.name]["position"] = Vector{}
		all_atoms[atom.name]["velocity"] = Vector{}
		all_atoms[atom.name]['potential'] = Vector{}
	}

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

		all_atoms[atoms[i].name]["position"] = atoms[i].position
		all_atoms[atoms[i].name]["velocity"] = atoms[i].velocity
		all_atoms[atoms[i].name]["potential"] = Vector.new(potential, 0.0, 0.0)
		
	}

	mut z := 1

	save_to_csv(all_atoms, state, 0)

	for t < simulaton_time {

		for i in 0..atoms.len {
			name := atoms[i].name
			mut pos := verlet_method_position(atoms[i], accelerations[i], dt)
			pos = validate_boundary(mut pos, size)

			atoms[i].position = pos

			all_atoms[name]["position"] = atoms[i].position
		}

		mut accelerations_new := []Vector{len: atoms.len}

		mut computed_values := map[string][]Vector{}
		
		for i in 0..atoms.len {
			name := atoms[i].name

			mut potential := Vector.new(0.0, 0.0, 0.0)
			mut new_force := Vector.new(0.0, 0.0, 0.0)

			for j in 0..atoms.len {
				if i == j {
					continue
				}

				mut force_ := Vector.new(0.0, 0.0, 0.0)
				mut potential_:= Vector.new(0.0, 0.0, 0.0)
				if '${j},${i}' in computed_values {
					c := computed_values['${j},${i}']
					potential_ = c[0]
					force_ = c[1].mul(-1)
				} else {
					potential__, force__ := lennard_jones(atoms[i].position, atoms[j].position, size)	
					potential = Vector.new(potential__, 0.0, 0.0)
					force_ = force__
					computed_values['${i},${j}'] = [potential_, force_]
				}

				potential += potential_ 
				new_force += force_
			}


			atoms[i].velocity = verlet_method_velocity(atoms[i], dt, new_force.div(atom_mass), accelerations[i])
			
			all_atoms[name]["potential"] = potential

			accelerations_new[i] = new_force.div(atom_mass)
		}

		mut velocities := []Vector{len: atoms.len}

		for i in 0..atoms.len {
			velocities[i] = atoms[i].velocity
		}

		if t < 1e-15 {
			temperature := kinetic_temperature(velocities)

			mut new_velocites := []Vector{len: atoms.len}

			for i in 0..atoms.len {
				new_velocites[i] = atoms[i].velocity.mul(math.sqrt(94.4 / temperature))
			}
			velocities = new_velocites.clone()
			new_velocites.clear()
		}


		for i in 0..atoms.len { 
			all_atoms[atoms[i].name]["velocity"] = velocities[i]
		}

		accelerations = accelerations_new.clone()
		accelerations_new.clear()

		println("Completed ${(t / simulaton_time * 100):.2}%")
		t += dt

		save_to_csv(all_atoms, state, z)
		z += 1	
	}

}

fn main() {
	run_simulation(State.lattice, 5)
}
