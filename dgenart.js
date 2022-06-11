function dgenart() {
  let _this = this;
  this.math = class dMath {
    static epsilon = 0.0001;
  
    static PI = 3.14159265359;
    static TAU = 6.28318530717;
    
    static lerp(a, b, t) {
      let d = b - a;
      return a + t*d;
    }
  
    static smooth_lerp(a, b, t) {
      let d = b - a;
      let t0 = t*t*t*(t*(t*6 - 15) + 10);
      return a + t0*d;
    }
  
    static eerp(a, b, t) {
      return Math.pow(a, t)*Math.pow(b, 1-t);
    }
  
    static loop(v, n) {
      let v0 = v;
      while (v0 < 0) { v0 += n; }
      while (v0 >= n) { v0 -= n; }
      return v0;
    }
  
    static loop2(v, n1, n2) {
      let v0 = v.copy();
      while (v0.x < 0) { v0.x += n1; }
      while (v0.x > n1) { v0.x -= n1; }
      while (v0.y < 0) { v0.y += n2; }
      while (v0.y > n2) { v0.y -= n2; }
      return v0;
    }
  
    static loop3(v, n1, n2, n3) {
      let v0 = v.copy();
      while (v0.x < 0) { v0.x += n1; }
      while (v0.x > n1) { v0.x -= n1; }
      while (v0.y < 0) { v0.y += n2; }
      while (v0.y > n2) { v0.y -= n2; }
      while (v0.z < 0) { v0.z += n3; }
      while (v0.z > n3) { v0.z -= n3; }
      return v0;
    }
  }
  this.vector2 = class dvector {
    constructor(x, y) {
      this.x = x;
      this.y = y;
    }
    static zero = new dvector(0, 0);
    static up = new dvector(0, 1);
    static down = new dvector(0, -1);
    static right = new dvector(1, 0);
    static left = new dvector(-1, 0);

    mag() { return Math.sqrt(this.x*this.x + this.y*this.y); }
    magsq() { return this.x*this.x + this.y*this.y; }
    heading() { return Math.atan2(this.y, this.x); }

    copy() { return new dvector(this.x, this.y); }

    normalize() { return this.div(this.mag()); }
    normalize_equals() { this.div_equals(this.mag()); }

    add(v) { return new dvector(this.x + v.x, this.y + v.y); }
    add_equals(v) { this.x += v.x; this.y += v.y; }

    sub(v) { return new dvector(this.x - v.x, this.y - v.y); }
    sub_equals(v) { this.x -= v.x; this.y -= v.y; }

    mult(n) { return new dvector(this.x * n, this.y * n); }
    mult_equals(n) { this.x *= n; this.y *= n; }

    mult_vec(v) { return new dvector(this.x * v.x, this.y * v.y); }
    mult_vec_equals(v) { this.x *= v.x; this.y *= v.y; }

    div(n) { return new dvector(this.x/n, this.y/n); }
    div_equals(n) { this.x /= n; this.y /= n; }

    div_vec(v) { return new dvector(this.x/v.x, this.y/v.y); }
    div_vec_equals(v) { this.x /= v.x; this.y /= v.y; }

    dot(v) { return(this.x*v.x + this.y*v.y) }

    rotate(a) {
      let x = Math.cos(a)*this.x - Math.sin(a)*this.y;
      let y = Math.sin(a)*this.x + Math.cos(a)*this.y;
      return new dvector(x, y);
    }
    rotate_equals(a) {
      let x = Math.cos(a)*this.x - Math.sin(a)*this.y;
      let y = Math.sin(a)*this.x + Math.cos(a)*this.y;
      this.x = x;
      this.y = y;
    }

    lerp(v, t) {
      let diff = v.sub(this);
      return this.add(diff.mult(t));
    }
    lerp_equals(v, t) {
      let diff = v.sub(this);
      this.add_equals(diff.mult(t));
    }

    reflect_x() {
      let v = this.copy();
      v.y *= -1;
      return v;
    }
    reflect_x_equals() {
      this.y *= -1;
    }
    reflect_y() {
      let v = this.copy();
      v.x *= -1;
      return v;
    }
    reflect_y_equals() {
      this.x *= -1;
    }
    reflect(a, b) {
      let pt = this.sub(a);
      let bt = b.sub(a);
      let pd = pt.dot(bt)/bt.mag();
      let pn = pt.sub(bt.normalize().mult(pd));
      let pr = pt.sub(pn.mult(2));
      let ptr = pr.add(a);
      return ptr;
    }
    reflect_equals(a, b) {
      let pt = this.sub(a);
      let bt = b.sub(a);
      let pd = pt.dot(bt)/bt.mag();
      let pn = pt.sub(bt.normalize().mult(pd));
      let pr = pt.sub(pn.mult(2));
      let ptr = pr.add(a);
      this.x = ptr.x; this.y = ptr.y
    }
  
    rot_lerp(v, t) {
      let m = this.mag();
      let n = this.lerp(v, t).normalize().mult(m);
      return n;
    }
    rot_lerp_equals(v, t) {
      let m = this.mag();
      this.lerp_equals(v, t);
      this.normalize_equals();
      this.mult_equals(m);
    }
  }
  this.perlin2 = class dPerlin2 {
    constructor(lookup_length, seed=dRandom.random(-1000, 1000)) {
      this.lookup = [];
      for (let i = 0; i < lookup_length; i++) { this.lookup.push(i); }
      this.rangen = dRandom.sfc32(seed, -seed/2, seed/5, -seed/10);
      dRandom.shuffle_seeded(this.lookup, this.rangen);
    }
  
    loop(n) {
      while (n >= this.lookup.length) {
        n -= this.lookup.length;
      }
      while (n < 0) {
        n += this.lookup.length;
      }
      return n;
    }
  
    hash(v) {
      return this.lookup[this.lookup[this.loop(v.x)] + this.loop(v.y)];
    }
  
    grad(n) {
      let v = n % 4;
      let x, y;
      if (v == 0 || v == 2) { x = -1; } else { x = 1; }
      if (v == 0 || v == 1) { y = -1; } else { y = 1; }
      return new dvector(x, y);
    }
  
    perlin_base(v) {
      let v0 = new dvector(Math.floor(v.x+dMath.epsilon), Math.floor(v.y+dMath.epsilon));
      let v1 = new dvector(Math.floor(v.x+dMath.epsilon), Math.ceil(v.y+dMath.epsilon));
      let v2 = new dvector(Math.ceil(v.x+dMath.epsilon), Math.floor(v.y+dMath.epsilon));
      let v3 = new dvector(Math.ceil(v.x+dMath.epsilon), Math.ceil(v.y+dMath.epsilon));
  
      let g0 = this.grad(this.hash(v0)).dot(v.sub(v0));
      let g1 = this.grad(this.hash(v1)).dot(v.sub(v1));
      let g2 = this.grad(this.hash(v2)).dot(v.sub(v2));
      let g3 = this.grad(this.hash(v3)).dot(v.sub(v3));
  
      let xt = v.sub(v0).x;
      let yt = v.sub(v0).y;
      let g01 = dMath.smooth_lerp(g0, g1, yt);
      let g23 = dMath.smooth_lerp(g2, g3, yt);
      let g0123 = dMath.smooth_lerp(g01, g23, xt);
      return g0123;
    }
  
    perlin(v, o=dvector.zero, f=1, a=1) {
      return a*this.perlin_base(v.mult(f/100).add(o));
    }
  
    perlin_octaves(v, p, n, o=dvector.zero, f=1, a=1) {
      let frequency = f;
      let amplitude = a;
      let total = 0;
      let max = 0;
      for (let i = 0; i < n; i++) {
        total += this.perlin(v, o, frequency, amplitude);
        max += amplitude;
        amplitude *= p;
        frequency *= 2;
      }
      return a*total/max;
    }
  
    perlin_warp(v, n, o=dvector.zero, f=1, a=1) {
      let t = 0;
      let values = [];
      for (let i = 0; i < n; i++) {
        t = this.perlin(v.add(new dvector(t, t), o, f, a)); 
        values.push(t);
      }
      return values;
    }
  
    perlin_warp_octaves(v, n, p=0.5, no=4, o=dvector.zero, f=1, a=1) {
      let t = 0;
      let values = [];
      for (let i = 0; i < n; i++) {
        t = this.perlin_octaves(v.add(new dvector(t, t), p, no, 0, f, a)); 
        values.push(t);
      }
      return values;
    }
  }
  this.random = class dRandom {
    constructor() {}
  
    static sfc32(a, b, c, d) {
      return function() {
        a >>>= 0; b >>>= 0; c >>>= 0; d >>>= 0; 
        var t = (a + b) | 0;
        a = b ^ b >>> 9;
        b = c + (c << 3) | 0;
        c = (c << 21 | c >>> 11);
        d = d + 1 | 0;
        t = t + d | 0;
        c = c + t | 0;
        return (t >>> 0) / 4294967296;
      }
    }
  
    static random(a=0, b=1) {
      let r = b-a;
      return a + Math.random()*r;
    }
  
    static random_seeded(rangen, a=0, b=1) {
      let r = b-a;
      return a + rangen()*r;
    }
  
    static shuffle(list) {
      let index = list.length;
      let ranindex;
      while (index != 0) {
        ranindex = Math.floor(this.random(0, index));
        index--;
        [list[index], list[ranindex]] = [list[ranindex], list[index]];
      }
      return list; 
    }
  
    static shuffle_seeded(list, rangen) {
      let index = list.length;
      let ranindex;
      while (index != 0) {
        ranindex = Math.floor(this.random_seeded(rangen, 0, index));
        index--;
  
        [list[index], list[ranindex]] = [list[ranindex], list[index]];
      }
      return list; 
    }
  
    static perlin_gen(seed=dRandom.random(-1000, 1000)) {
      return new dPerlin2(256, seed);
    }
  }
}