module vector
import math

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

pub fn Vector.zero() Vector {
	return Vector { x: 0.0, y: 0.0, z: 0.0 }
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