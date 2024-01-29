import math

type Vector = [3]f64

struct Atom {
pub:
	mass: f64 = 39.948 * 931.394 * math.power(10, 9) / (299_792_458 * 299_792_458)
pub mut:
	position: Vector
	velocity: Vector
}

fn (v Vector) + (w Vector) Vector {
	return [v[0] + w[0], v[1] + w[1], v[2] + w[2]]
}

fn (v Vector) - (w Vector) Vector {
	return [v[0] - w[0], v[1] - w[1], v[2] - w[2]]
}

fn (v Vector) norm() f64 {
	return math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
}

fn (v Vector) scale(scalar f64) Vector {
	return [v[0]*scalar, v[1]*scalar, v[2]*scalar]
}

fn (v Vector) dot(w Vector) f64 {
	return v[0]*w[0] + v[1]*w[1] + v[2]*w[2]
}

fn VerletPosition(mut atom Atom, dt f64, force Vector) {
	atom.position = atom.position + atom.velocity.scale(dt) + force.scale(0.5*dt*dt*atom.mass)
}

fn VerletVelocity(mut atom Atom, dt f64, newforce Vector, oldForce Vector) {
	atom.velocity = atom.velocity + (newforce + oldForce).scale(0.5*dt/atom.mass) 
}

fn initalize() []Atom {
	return [
		Atom{
			Vector(0.0, 0.0, 0.0),
			Vector(0.0, 0.0, 0.0),
		},
		Atom {
			Vector(4.0, 0.0, 0.0),
			Vector(0.0, 0.0, 0.0),
		}
	]
}

fn main() {
	atoms = initalize()
}


