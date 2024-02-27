module vector
import math

// Vector struct, for 3D vectors
pub struct Vector {
pub mut:
	x f64
	y f64
	z f64
}

// Create a new vector
// : x {f64} x component
// : y {f64} y component
// : z {f64} z component
// ; {Vector} the new vector
pub fn Vector.new(x f64, y f64, z f64) Vector {
	return Vector { x: x, y: y, z: z }
}

// Create a new vector zero vector
// ; {Vector} the new vector
pub fn Vector.zero() Vector {
	return Vector { x: 0.0, y: 0.0, z: 0.0 }
}

@[inline]
// Add two vectors
// : v {Vector} the first vector
// : u {Vector} the second vector
// ; {Vector} the sum of the two vectors
pub fn (v Vector) + (u Vector) Vector {
	return Vector { x: v.x + u.x, y: v.y + u.y, z: v.z + u.z }
}


@[inline]
// Subtract two vectors
// : v {Vector} the first vector
// : u {Vector} the second vector
// ; {Vector} the difference of the two vectors
pub fn (v Vector) - (u Vector) Vector {
	return Vector { x: v.x - u.x, y: v.y - u.y, z: v.z - u.z }
}

// Multiply a vector by a scalar
// : v {Vector} the vector
// : s {f64} the scalar
// ; {Vector} the product of the vector and the scalar
pub fn (v Vector) mul(s f64) Vector {
	return Vector { x: v.x * s, y: v.y * s, z: v.z * s }
}

// Computes the norm of a vector
// : v {Vector} the vector
// ; {f64} the norm of the vector
pub fn (v Vector) norm() f64 {
	return math.sqrt(v.x * v.x + v.y * v.y + v.z * v.z)
}

// Divides a vector by a scalar
// : v {Vector} the vector
// : s {f64} the scalar
// ; {Vector} the division of the vector by the scalar
pub fn (v Vector) div(s f64) Vector {
	if s == 0.0 {
		panic("Division by zero")
	}

	return Vector { x: v.x / s, y: v.y / s, z: v.z / s }
}

// Get the i:th component of a vector
// : v {Vector} the vector
// : i {int} the index
// ; {f64} the i:th component of the vector
pub fn (v Vector) get(i int) ?f64 {
	return match i {
		0 { v.x }
		1 { v.y }
		2 { v.z }
		else { none }
	}
}

// Get vector's components as an array
// : v {Vector} the vector
// ; {[]f64} the vector's components
pub fn (v Vector) to_array() []f64 {
	return [v.x, v.y, v.z]
}