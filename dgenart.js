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
  this.vector3 = class dvector3 {
    constructor(x, y, z) {
      this.x = x;
      this.y = y;
      this.z = z;
    }
    static zero = new dvector3(0, 0, 0);
    static up = new dvector3(0, 1, 0);
    static down = new dvector3(0, -1, 0);
    static left = new dvector3(-1, 0, 0);
    static right = new dvector3(1, 0, 0);
    static forward = new dvector3(0, 0, 1);
    static backward = new dvector3(0, 0, -1);
    mag() {
      return Math.sqrt(this.x*this.x + this.y*this.y + this.z*this.z);
    }
    magsq() {
      return this.x*this.x + this.y*this.y + this.z*this.z;
    }
    copy() {
      return new dvector3(this.x, this.y, this.z);
    }
    normalize() {
      return this.div(this.mag());
    }
    normalize_equals() {
      this.div_equals(this.mag());
    }
    add(v) {
      return new dvector3(this.x + v.x, this.y + v.y, this.z + v.z);
    }
    add_equals(v) {
      this.x += v.x;
      this.y += v.y;
      this.z += v.z;
    }
    sub(v) {
      return new dvector3(this.x - v.x, this.y - v.y, this.z - v.z);
    }
    sub_equals(v) {
      this.x -= v.x;
      this.y -= v.y;
      this.z -= v.z;
    }
    mult(n) {
      return new dvector3(this.x * n, this.y * n, this.z * n);
    }
    mult_equals(n) {
      this.x *= n;
      this.y *= n;
      this.z *= n;
    }
    mult_vec(v) {
      return new dvector3(this.x * v.x, this.y * v.y, this.z * v.z);
    }
    mult_vec_equals(v) {
      this.x *= v.x;
      this.y *= v.y;
      this.z *= v.z;
    }
    div(n) {
      return new dvector3(this.x/n, this.y/n, this.z/n);
    }
    div_equals(n) {
      this.x /= n;
      this.y /= n;
      this.z /= n;
    }
    div_vec(v) {
      return new dvector3(this.x/v.x, this.y/v.y, this.z/v.z);
    }
    div_vec_equals(v) {
      this.x /= v.x;
      this.y /= v.y;
      this.z /= v.z;
    }
    lerp(v, t) {
      let diff = v.sub(this);
      return this.add(diff.mult(t));
    }
    lerp_equals(v, t) {
      let diff = v.sub(this);
      this.add_equals(diff.mult(t));
    }
    dot(v) {
      return(this.x*v.x + this.y*v.y + this.z*v.z);
    }
    cross(v) {
      let x = this.y*v.z - this.z-v.y;
      let y = this.z*v.x - this.x-v.z;
      let z = this.x*v.y - this.y-v.x;
      return new dvector3(x, y, z);
    }
    rotate_axis_angle(axis, angle) {
      let v = axis.normalize();
      let qx = v.x*Math.sin(angle/2);
      let qy = v.y*Math.sin(angle/2);
      let qz = v.z*Math.sin(angle/2);
      let qw = Math.cos(angle/2);
      let q = new dquaternion(qx, qy, qz, qw).normalize();
      console.log(q);
      console.log(q.to_euler());
      return this.rotate_quaternion(q);
    }
    rotate_quaternion(q) {
      let v = this.copy();
  
      var tx = 2 * (q.y * v.z - q.z * v.y);
      var ty = 2 * (q.z * v.x - q.x * v.z);
      var tz = 2 * (q.x * v.y - q.y * v.x);
  
      return new dvector3(
        v.x + q.w * tx + q.y * tz - q.z * ty,
        v.y + q.w * ty + q.z * tx - q.x * tz,
        v.z + q.w * tz + q.x * ty - q.y * tx);
    }
  }
  this.geom = class dGeom {
    constructor() {}
  
    static av_points(p) {
      let av = dvector.zero;
      for (let i = 0; i < p.length; i++) {
        av.add_equals(p[i]);
      }
      return av.div(p.length);
    }
    
    static cubic_lerp(a, b, c, d, t) {
      let ab = a.lerp(b, t);
      let bc = b.lerp(c, t);
      let cd = c.lerp(d, t);
      let abc = ab.lerp(bc, t);
      let bcd = bc.lerp(cd, t);
      return abc.lerp(bcd, t);
    }
  
    static orientation(a, b, c) {
      let val = (b.y - a.y) * (c.x - b.x) -
                (b.x - a.x) * (c.y - b.y);
      if (val == 0) { return 0; } //colinear
      else if (val > 0) { return 1; } //clockwise
      else { return 2; } //counterclockwise
    }
  
    static on_segment(p, q1, q2) {
      if (p.x <= Math.max(q1.x, q2.x) &&
          p.x >= Math.min(q1.x, q2.x) &&
          p.y >= Math.max(q1.y, q2.y) &&
          p.y <= Math.min(q1.y, q2.y)) { return true; }
      else { return false; }
    }
  
    static line_intersect(p1, p2, q1, q2) {
      let o1 = this.orientation(p1, p2, q1);
      let o2 = this.orientation(p1, p2, q2);
      let o3 = this.orientation(q1, q2, p1);
      let o4 = this.orientation(q1, q2, p2);
  
      if (o1 != o2 && o3 != o4) { return true; }
      if (o1 == 0 && this.on_segment(q1, p1, p2)) { return true; }
      if (o2 == 0 && this.on_segment(q2, p1, p2)) { return true; }
      if (o3 == 0 && this.on_segment(p1, q1, q2)) { return true; }
      if (o4 == 0 && this.on_segment(p2, q1, q2)) { return true; }
  
      return false;
    }
  
    static point_in_poly(point, poly, transform=dvector.zero, scale=1) {
      let bounds = poly.get_bounds();
      let points = poly.get_points(transform, scale)
      let p1 = point;
      let p2 = new dvector(bounds[1].x*10*scale, point.y);
      let num_intersect = 0;
      for (let i = 0; i < points.length; i++) {
        let q1 = points[i];
        let q2i = i + 1; if (q2i >= points.length) { q2i = 0; }
        let q2 = points[q2i];
        if (this.line_intersect(p1, p2, q1, q2)) { num_intersect++; }
      }
      if (num_intersect % 2 == 0) { return false; }
      else { return true; }
    }
  
    static dist_to_plane(point, plane_point, plane_normal) {
      let v = point.sub(plane_point);
      let dist = v.dot(plane_normal);
      return dist;
    }
  }
  this.phys2 = class dPhys {
    constructor() {}
  
    static grav(a, b, g, dt=1, affect=2, scale=1) {
      let vec = b.pos.sub(a.pos);
      if (vec.mag() < 1) { return; }
      let f = dt*(g * a.mass * b.mass)/(vec.magsq());
      if (affect == 0 || affect == 2) {
        a.add_force(vec.normalize().mult(f));
      }
      if (affect == 1 || affect == 2) {
        b.add_force(vec.normalize().mult(-1*f));
      }
      return vec.normalize().mult(f);
    }
  }
  this.grid2 = class dGrid {
    constructor(w, h, init=0, loop=false) {
      this.width = w;
      this.height = h;
      this.init = init;
      this.loop = loop;
      this.grid = [];
      for (let i = 0; i < w; i++) {
        this.grid.push([]);
        for (let j = 0; j < h; j++) {
          this.grid[i].push(init);
        }
      }
    }
  
    get(x, y, b=null) {
      if (x < 0 || y < 0 || x >= this.width || y >= this.height) {
        if (this.loop) {
          let tx = x; let ty = y;
          while (tx < 0) { tx += this.width; }
          while (ty < 0) { ty += this.height; }
          while (tx >= this.width) { tx -= this.width; }
          while (ty >= this.height) { ty -= this.height; }
          return this.grid[tx][ty];
        }
        else if (b == null) { return this.init; }
        else { return b; }
      }
      else if (x != NaN && y != NaN) { return this.grid[x][y]; }
      else { return this.grid[0][0]; }
    }
  
    set(x, y, v) {
      if (x >= 0 && y >= 0 && x < this.width && y < this.height) {
        this.grid[x][y] = v;
      }
      else if (this.loop) {
        let tx = x; let ty = y;
          while (tx < 0) { tx += this.width; }
          while (ty < 0) { ty += this.height; }
          while (tx >= this.width) { tx -= this.width; }
          while (ty >= this.height) { ty -= this.height; }
          this.grid[tx][ty] = v;
      }
    }
  
    interpolate(v, o=dvector.zero, w=this.width, h=this.height) {
      let vt = v.sub(o);
      let vf = new dvector(Math.floor(this.width*vt.x/w), Math.floor(this.height*vt.y/h));
      let v0 = vt.sub(vf);
      let v1 = vt.sub(vf.add(new dvector(1, 0)));
      let v2 = vt.sub(vf.add(new dvector(0, 1)));
      let v3 = vt.sub(vf.add(new dvector(1, 1)));
      let va = v0.mag() + v1.mag() + v2.mag() + v3.mag();
      let r0 = [this.get(vf.x, vf.y), v0.mag()/va];
      let r1 = [this.get(vf.x + 1, vf.y), v1.mag()/va];
      let r2 = [this.get(vf.x, vf.y + 1), v2.mag()/va];
      let r3 = [this.get(vf.x + 1, vf.y + 1), v3.mag()/va];
      return [r0, r1, r2, r3];
    }
  }
  this.grid3 = class dGrid3 {
    constructor(w, h, d, init=0, loop=false) {
      this.width = w;
      this.height = h;
      this.depth = d;
      this.init = init;
      this.loop = loop;
      this.grid = [];
      for (let i = 0; i < w; i++) {
        this.grid.push([]);
        for (let j = 0; j < h; j++) {
          this.grid[i].push([]);
          for (let k = 0; k < d; k++) {
            this.grid[i][j].push(init);
          }
        }
      }
    }
  
    get(x, y, z, b=null) {
      if (x < 0 || y < 0 || z < 0 || x >= this.width || y >= this.height || z >= this.depth) {
        if (this.loop) {
          let tx = x; let ty = y; let tz = z;
          while (tx < 0) { tx += this.width; }
          while (ty < 0) { ty += this.height; }
          while (tz < 0) { tz += this.depth; }
          while (tx >= this.width) { tx -= this.width; }
          while (ty >= this.height) { ty -= this.height; }
          while (tz >= this.height) { tz -= this.depth; }
          return this.grid[tx][ty][tz];
        }
        else if (b == null) { return this.init; }
        else { return b; }
      }
      else { return this.grid[x][y][z]; }
    }
  
    set(x, y, z, v) {
      if (x >= 0 && y >= 0 && z >= 0 && x < this.width && y < this.height && z < this.depth) {
        this.grid[x][y][z] = v;
      }
      else if (this.loop) {
        let tx = x; let ty = y; let tz = z;
        while (tx < 0) { tx += this.width; }
        while (ty < 0) { ty += this.height; }
        while (tz < 0) { tz += this.depth; }
        while (tx >= this.width) { tx -= this.width; }
        while (ty >= this.height) { ty -= this.height; }
        while (tz >= this.height) { tz -= this.depth; }
        this.grid[tx][ty][tz] = v;
      }
    }
  }
  this.grid2_buffer = class dGrid_Buffer {
    constructor(width, height, init=0, loop=false) {
      this.width = width;
      this.height = height;
      this.init = init;
      this.loop = loop;
      this.grid = new dGrid(width, height, init, loop);
      this.grid0 = new dGrid(width, height, init, loop);
    }
    get(x, y, b=null) {
      return this.grid.get(x, y, b);
    }
    get_both(x, y, b=null) {
      return [this.grid.get(x, y, b), this.grid0.get(x, y, b)];
    }
    set(x, y, v) {
      this.grid0.set(x, y, v);
    }
    set_both(x, y, v) {
      this.grid.set(x, y, v);
      this.grid0.set(x, y, v);
    }
    buffer() {
      for (let i = 0; i < this.width; i++) {
        for (let j = 0; j < this.height; j++) {
          this.grid.set(i, j, this.grid0.get(i, j));
        }
      }
    }
    reverse_buffer() {
      for (let i = 0; i < this.width; i++) {
        for (let j = 0; j < this.height; j++) {
          this.grid0.set(i, j, this.grid.get(i, j));
        }
      }
    }
  }
  this.grid3_buffer = class dGrid3_Buffer {
    constructor(width, height, depth, init=0, loop=false) {
      this.width = width;
      this.height = height;
      this.depth = depth;
      this.init = init;
      this.loop = loop;
      this.grid = new dGrid3(width, height, depth, init, loop);
      this.grid0 = new dGrid3(width, height, depth, init, loop);
    }
    get(x, y, z, b=null) {
      return this.grid.get(x, y, z, b);
    }
    get_both(x, y, z, b=null) {
      return [this.grid.get(x, y, z, b), this.grid0.get(x, y, z, b)];
    }
    set(x, y, z, v) {
      this.grid0.set(x, y, z, v);
    }
    set_both(x, y, z, v) {
      this.grid.set(x, y, z, v);
      this.grid0.set(x, y, z, v);
    }
    buffer() {
      for (let x = 0; x < this.width; x++) {
        for (let y = 0; y < this.height; y++) {
          for (let z = 0; z < this.depth; z++) {
            this.grid.set(x, y, z, this.grid0.get(x, y, z));
          }
        }
      }
    }
    reverse_buffer() {
      for (let x = 0; x < this.width; x++) {
        for (let y = 0; y < this.height;yj++) {
          for (let z = 0; z < this.depth; z++) {
            this.grid0.set(x, y, z, this.grid.get(x, y, z));
          }
        }
      }
    }
  }
  this.quaternion = class dquaternion {
    constructor(x, y, z, w) {
      this.x = x;
      this.y = y;
      this.z = z;
      this.w = w;
    }
    static zero = new dquaternion(0, 0, 0, 0);
  
    mag() {
      return Math.sqrt(this.x*this.x + this.y*this.y + this.z*this.z + this.w*this.w);
    }
    magsq() {
      return this.x*this.x + this.y*this.y + this.z*this.z  + this.w*this.w;
    }
    copy() {
      return new dquaternion(this.x, this.y, this.z, this.w);
    }
    normalize() {
      return this.div(this.mag());
    }
    normalize_equals() {
      this.div_equals(this.mag());
    }
    add(v) {
      return new dquaternion(this.x + v.x, this.y + v.y, this.z + v.z, this.w + v.w);
    }
    add_equals(v) {
      this.x += v.x;
      this.y += v.y;
      this.z += v.z;
      this.w += v.w;
    }
    sub(v) {
      return new dquaternion(this.x - v.x, this.y - v.y, this.z - v.z, this.w - v.w);
    }
    sub_equals(v) {
      this.x -= v.x;
      this.y -= v.y;
      this.z -= v.z;
      this.w -= v.w;
    }
    mult(n) {
      return new dquaternion(this.x * n, this.y * n, this.z * n, this.w * n);
    }
    mult_equals(n) {
      this.x *= n;
      this.y *= n;
      this.z *= n;
      this.w *= n;
    }
    mult_vec(v) {
      return new dquaternion(this.x * v.x, this.y * v.y, this.z * v.z, this.w * v.w);
    }
    mult_vec_equals(v) {
      this.x *= v.x;
      this.y *= v.y;
      this.z *= v.z;
      this.w *= v.w;
    }
    div(n) {
      return new dquaternion(this.x/n, this.y/n, this.z/n, this.w/n);
    }
    div_equals(n) {
      this.x /= n;
      this.y /= n;
      this.z /= n;
      this.w /= n;
    }
    div_vec(v) {
      return new dquaternion(this.x/v.x, this.y/v.y, this.z/v.z, this.w/v.w);
    }
    div_vec_equals(v) {
      this.x /= v.x;
      this.y /= v.y;
      this.z /= v.z;
      this.w /= v.w;
    }
    to_euler() {
      let sinr_cosp = 2 * (this.w * this.x + this.y *this.z);
      let cosr_cosp = 1 - 2 * (this.x * this.x + this.y * this.y);
      let x = Math.atan2(sinr_cosp, cosr_cosp);
  
      let sinp = 2 * (this.w * this.y - this.z * this.x);
      let y;
      if (Math.abs(sinp) >= 1) {
        y = Math.PI/2 * Math.sign(sinp);
      }
      else {
        y = Math.asin(sinp);
      }
  
      let siny_cosp = 2 * (this.w * this.z + this.x * this.y);
      let cosy_cosp = 1 - 2 * (this.y * this.y * this.z * this.z);
      let z = Math.atan2(siny_cosp, cosy_cosp);
  
      return new dvector3(x, y, z);
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
  this.camera = class dCamera {
    constructor(pos, dir, up, length, width, height) {
      this.pos = pos;
      this.dir = dir;
      this.up = up;
      this.length = length;
      this.width = width;
      this.height = height;
      this.corners = [];
      this.screen_centre = this.pos.add(this.dir.normalize().mult(this.length));
      this.right = this.dir.cross(this.up).normalize();
      this.corners[0] = this.screen_centre.add(this.right.mult(-width/2)).add(this.up.mult(height/2));
      this.corners[1] = this.screen_centre.add(this.right.mult(width/2)).add(this.up.mult(height/2));
      this.corners[2] = this.screen_centre.add(this.right.mult(-width/2)).add(this.up.mult(-height/2));
      this.corners[3] = this.screen_centre.add(this.right.mult(width/2)).add(this.up.mult(-height/2));
    }
  
    regen() {
      this.corners = [];
      this.screen_centre = this.pos.add(this.dir.normalize().mult(this.length));
      this.right = this.dir.cross(this.up).normalize();
      this.corners[0] = this.screen_centre.add(this.right.mult(-width/2)).add(this.up.mult(height/2));
      this.corners[1] = this.screen_centre.add(this.right.mult(width/2)).add(this.up.mult(height/2));
      this.corners[2] = this.screen_centre.add(this.right.mult(-width/2)).add(this.up.mult(-height/2));
      this.corners[3] = this.screen_centre.add(this.right.mult(width/2)).add(this.up.mult(-height/2));
    }
  
    get_xy(point) {
      let x = dGeom.dist_to_plane(point, this.corners[2], this.right);
      let y = dGeom.dist_to_plane(point, this.corners[2], this.up);
      return new dvector(x/this.width, y/this.height);
    }
  
    get_isect(point) {
      let u = point.sub(this.pos);
      let dot = this.dir.dot(u);
  
      if (Math.abs(dot) > 0.0001) {
        let w = this.pos.sub(this.screen_centre);
        let fac = this.dir.dot(w)/dot;
        u.mult(fac);
        return this.pos.add(u);
      }
      else {
        return this.screen_centre;
      }
    }
  
    get_screen_coords(point) {
      let coords = this.get_isect(point);
      let dist = dGeom.dist_to_plane(point, this.screen_centre, this.dir);
      let xy = this.get_xy(coords);
      return new dvector3(xy.x, xy.y, dist/this.length);
    }
  
  }
  this.poly = class dPoly {
    constructor(points) {
      this.points = [...points];
      for (let i = 0; i < this.points.length; i++) {
        this.points[i] = this.points[i].copy();
      }
    }
  
    static unit_triangle = new dPoly([
      new dvector(1, 0),
      new dvector(1, 0).rotate(2*Math.PI/3),
      new dvector(1, 0).rotate(-2*Math.PI/3)
    ]);
  
    static unit_square = new dPoly([
      new dvector(1, 0).rotate(Math.PI/4),
      new dvector(1, 0).rotate(Math.PI/2 + Math.PI/4),
      new dvector(1, 0).rotate(Math.PI + Math.PI/4),
      new dvector(1, 0).rotate(-Math.PI/2 + Math.PI/4)
    ]);
  
    static circle(r, n) {
      let points = [];
      let v = new dvector(r, 0);
      for (let i = 0; i < n; i++) {
        points.push(v.rotate(2*Math.PI*(i/n)));
      }
      return new dPoly(points);
    }
  
    copy() {
      return new dPoly([...this.points]);
    }
  
    get_centre() {
      let av = new dvector(0, 0);
      for (let i = 0; i < this.points.length; i++) {
        av.add_equals(this.points[i]);
      }
      av.div_equals(this.points.length);
      return av;
    }
  
    get_bounds() {
      let minx = [0, 10000];
      let miny = [0, 10000];
      let maxx = [0, -10000];
      let maxy = [0, -10000];
      for (let i = 0; i < this.points.length; i++) {
        let p = this.points[i];
        //console.log(p)
        if (p.x < minx[1]) { minx = [i, p.x]; }
        if (p.y < miny[1]) { miny = [i, p.y]; }
        if (p.x > maxx[1]) { maxx = [i, p.x]; }
        if (p.y > maxy[1]) { maxy = [i, p.y]; }
      }
      let minbounds = new dvector(minx[1], miny[1]);
      let maxbounds = new dvector(maxx[1], maxy[1]);
      return [minbounds, maxbounds];
    }
  
    get_points(transform=dvector.zero, scale=1) {
      let points = [];
      let centre = this.get_centre();
      for (let i = 0; i < this.points.length; i++) {
        let p = this.points[i].copy();
        let v = p.sub(centre);
        v.mult_equals(scale);
        p = centre.add(v);
        p.add_equals(transform);
        points.push(p);
      }
      return points;
    }
  
    move(v) {
      for (let i = 0; i < this.points.length; i++) {
        this.points[i].add_equals(v);
      }
    }
  
    scale(n) {
      let c = this.get_centre();
      for (let i = 0; i < this.points.length; i++) {
        let v = this.points[i].sub(c);
        v.mult_equals(n);
        this.points[i] = c.add(v);
      }
    }
  
    rotate(a) {
      let c = this.get_centre();
      for (let i = 0; i < this.points.length; i++) {
        let v = this.points[i].sub(c);
        v.rotate_equals(a);
        this.points[i] = c.add(v);
      }
    }
  
    distort(perlin, freq, amp, rx=false, ry=false) {
      let centre = this.get_centre();
      for (let i = 0; i < this.points.length; i++) {
        let p = this.points[i].copy();
        let v = p.sub(centre);
        let vp = v.copy();
        if (rx) { vp.y = Math.abs(vp.y); }
        if (ry) { vp.x = Math.abs(vp.x); }
        v.mult_equals(1 + perlin.perlin(vp, dvector.zero, freq, amp));
        this.points[i] = centre.add(v);
      }
    }
  
    get_points_distort(perlin, freq=1, amp=1, transform=dvector.zero, scale=1) {
      let centre = this.get_centre();
      let points = [];
      for (let i = 0; i < this.points.length; i++) {
        let p = this.points[i].copy();
        let v = p.sub(centre);
        v.mult_equals(1 + perlin.perlin(v, dvector.zero, freq, amp));
        v.mult_equals(scale)
        points.push(centre.add(v).add(transform));
      }
      return points; 
    }
  
    crinkle(n, m) {
      for (let i = 0; i < n; i++) {
        let current = this.get_points();
        let next = [];
        for (let j = 0; j < current.length; j++) {
          let jn = dMath.loop(j + 1, current.length);
          let vec = current[jn].sub(current[j]).div(2);
          vec.rotate_equals(Math.random()*m*2*Math.PI - m*Math.PI);
          next.push(current[j].add(vec));
        }
        let ret = [];
        for (let j = 0; j < current.length; j++) {
          ret.push(current[j]);
          ret.push(next[j]);
        }
        this.points = ret;
      }
    }
  
    nearest_point(p) {
      let min = [0, 1000000];
      for (let i = 0; i < this.points.length; i++) {
        let vec = this.points[i].sub(p);
        if (vec.mag() < min[1]) { min[0] = i; min[1] = vec.mag(); }
      }
      return(this.points[min[0]]);
    }
  
    fracture(pa, pb, ca=0.5, cb=0.5) {
      if (pb < pa) { let temp = pa; pa = pb; pb = temp; }
      let pa0 = dMath.loop(pa + 1, this.points.length);
      let pb0 = dMath.loop(pb + 1, this.points.length);
      let a = this.points[pa].lerp(this.points[pa0], ca);
      let b = this.points[pb].lerp(this.points[pb0], cb);
      let s1 = []; let s2 = [a];
      for (let i = 0; i <= pa; i++) {
        s1.push(this.points[i].copy());
      }
      s1.push(a); s1.push(b);
      for (let i = pa0; i <= pb; i++) {
        s2.push(this.points[i].copy());
      }
      s2.push(b);
      for (let i = pb0; i < this.points.length; i++) {
        s1.push(this.points[i].copy());
      }
      return([new dPoly(s1), new dPoly(s2)]);
    }
  
    proportional_move(n, v, m) {
      let o = this.points[n].copy();
      for (let i = 0; i < this.points.length; i++) {
        let p = this.points[i];
        if (p.sub(o).mag() < m) {
          if (i == n) { p.add_equals(v); }
          else { p.add_equals(v.mult( (m - p.sub(o).mag())/m )); }
        }
      }
    }
  }
  this.spline_point = class dSpline_Point {
    constructor(pos, fa, fm, ba, bm) {
      this.pos = pos;
      this.fa = fa;
      this.fm = fm;
      this.ba = ba;
      this.bm = bm;
    }
  
    get_forward() {
      let vec = new dvector(this.fm, 0);
      vec.rotate_equals(this.fa);
      return vec.add(this.pos);
    }
  
    get_back() {
      let vec = new dvector(this.bm, 0)
      vec.rotate_equals(this.ba);
      return vec.add(this.pos);
    }
  
    set_curved(normal, mag=1) {
      let dir = normal.rotate(dMath.PI/2);
      dir.normalize_equals().mult_equals(mag);
      this.fa = dir.heading();
      this.ba = dir.heading() + PI;
      this.fm = mag;
      this.bm = mag;
    }
  }
  this.spline = class dSpline {
    constructor(points) {
      this.points = points;
    }
  
    static from_points(points) {
      let spoints = [];
      for (let i = 0; i < points.length; i++) {
        let s = new dSpline_Point(points[i], 0, 1, 0, 1);
        spoints.push(s);
      }
      return new dSpline(spoints);
    }
  
    get_points(n, vec=null, loop=false) {
      let o;
      let points = [];
      if (vec == null) { o = dvector.zero; }
      else { o = vec.copy(); }
  
      let num_points = this.points.length - 1;
      if (loop) { num_points = this.points.length; }
      for (let p = 0; p < num_points; p++) {
        let pn = p + 1;
        if (pn >= this.points.length) { pn = 0; }
        let a = this.points[p].pos;
        let b = this.points[p].get_forward();
        let d = this.points[pn].pos;
        let c = this.points[pn].get_back();
  
        for (let i = 0; i < n; i++) {
          let t = i/n;
          points.push(dGeom.cubic_lerp(a, b, c, d, t).add(o));
        }
      }
  
      return points;
    }
  
    get_poly(n, vec=null, loop=false) {
      let p = this.get_points(n, vec, loop);
      return new dPoly(p);
    }
  }
  this.mover2 = class dMover {
    constructor(pos, vel, mass=1) {
      this.pos = pos;
      this.vel = vel;
      this.mass = mass;
    }
  
    update(dt=1, fr=0) {
      this.vel.mult_equals((1-(fr*dt)));
      this.pos.add_equals(this.vel.mult(dt));
    }
  
    add_force(f) {
      this.vel.add_equals(f.div(this.mass));
    }
  
    collide(collider, repel=3) {
      let collision = collider.collide(this.pos);
      if (!collision[0]) { return; }
      else { 
        this.pos = collision[1].add(collision[2].mult(repel)); 
        //this.pos.sub_equals(this.vel);
        this.vel.reflect_equals(dvector.zero, collision[2]); 
        this.vel.mult_equals(-1); 
        this.pos.add_equals(this.vel);
      }
    }
  }
  this.collider2 = class dCollider {
    constructor() {}
    
    static circle = class dCircleCollider {
      constructor(pos, radius) {
        this.pos = pos;
        this.radius = radius;
      }
  
      collide(p) {
        if (p.sub(this.pos).mag() <= this.radius) {
          let vec = p.sub(this.pos);
          return [true, vec.normalize().mult(this.radius), vec.normalize()];
        }
        else { return [false, dvector.zero, dvector.right]; }
      }
    }
  
    static rect = class dRectCollider {
      constructor(pos, side_length_a, side_length_b) {
        this.pos = pos;
        this.side_length_a = side_length_a;
        this.side_length_b = side_length_b;
      }
  
      corners() {
        let a = new dvector(-this.side_length_a/2, -this.side_length_b/2);
        let b = new dvector(this.side_length_a/2, -this.side_length_b/2);
        let c = new dvector(-this.side_length_a/2, this.side_length_b/2);
        let d = new dvector(this.side_length_a/2, this.side_length_b/2);
        return [a, b, c, d];
      }
  
      collide(p) {
        let corners = this.corners();
        if (p.x > corners[0].x && p.x < corners[3].x && p.y > corners[0].y && p.y < corners[3].y) {
          intersects = [];
          let right_dist = corners[3].x - p.x;
          let left_dist = p.x - corners[0].x;
          let up_dist = corners[3].y - p.y;
          let down_dist = p.y - corners[0].y;
          let min = Math.min(right_dist, Math.min(left_dist, Math.min(up_dist, down_dist)));
          if (min == right_dist) {
            return [true, new dvector(corners[3].x, p.y), dvector.right];
          }
          if (min == left_dist) {
            return [true, new dvector(corners[0].x, p.y), dvector.left];
          }
          if (min == up_dist) {
            return [true, new dvector(p.x, corners[3].y), dvector.up];
          }
          if (min == down_dist) {
            return [true, new dvector(p.x, corners[0].y), dvector.down];
          }
        } else { return [false, dvector.zero]; }
      }
    }
  
    static poly = class dPolyCollider {
      constructor(poly, pos=dvector.zero, scale=1) {
        this.poly = poly;
        this.pos = pos;
        this.scale = scale;
      }
  
      collide(p) {
        let c = dGeom.point_in_poly(p, this.poly, this.pos, this.scale);
        if (!c) { return [false, dvector.zero, dvector.right]; }
        else {
          let min = [0, 10000000000];
          let min0;
          let points = this.poly.get_points(this.pos, this.scale);
          for (let i = 0; i < points.length; i++) {
            let vec = points[i].sub(p);
            if (vec.mag() < min[1]) { min[0] = i; min[1] = vec.mag(); }
          }
          let minp = dMath.loop(min[0] + 1, points.length);
          let minn = dMath.loop(min[0] - 1, points.length);
          let card = -1;
          if (points[minp].sub(points[min[0]]).mag() < points[minn].sub(points[min[0]]).mag()) {
            min0 = minp;
          }
          else { min0 = minn; card = 1; }
          let bounds = this.poly.get_bounds();
          let bvec = bounds[1].sub(bounds[0]);
          let v = points[min0].sub(points[min[0]]);
          let u = points[min[0]].sub(p);
          let t = -1 * v.dot(u)/v.dot(v);
          if (t >= 1) { return [true, points[min0], points[min0].sub(this.poly.get_centre()).normalize()]; }
          else if (t <= 0) { return [true, points[min[0]], points[min[0]].sub(this.poly.get_centre()).normalize()]; }
          let p0 = points[min[0]].lerp(points[min0], t)
          let n = v.rotate(card*Math.PI/2).normalize();
          return[true, p0, n];
        }
      } 
    }
  }
  this.space_growth_node = class dSpace_Growth_Node {
    constructor(pos, parent=null) {
      this.pos = pos;
      this.parent = parent;
      this.children = [];
    }
  }
  this.space_growth_attractor = class dSpace_Growth_Attractor {
    constructor(pos, kill_radius) {
      this.pos = pos;
      this.kill_radius = kill_radius;
    }
  }
  this.space_growth_tree = class dSpace_Growth_Tree {
    constructor(origin, interval) {
      this.origin = origin;
      this.interval = interval;
      this.nodes = [];
      this.attractors = [];
      this.nodes.push(new dSpace_Growth_Node(origin.copy()));
    }
  
    attractor_gen_circle(r, n, kr) {
      for (let i = 0; i < n; i++) {
        let vec = new dvector(dRandom.random(-r, r), dRandom.random(-r, r));
        if (vec.mag() > r) { vec = vec.normalize().mult(r); }
        this.attractors.push(new dSpace_Growth_Attractor(vec.add(this.origin), kr));
      }
    }
  
    attractor_gen_spline(spline, n, kr, mult=1) {
      let points = spline.get_points(n*10, dvector.zero, true);
      for (let i = 0; i < n; i++) {
        let j = Math.floor(dRandom.random(0, points.length-1));
        this.attractors.push(new dSpace_Growth_Attractor(points[j].mult(mult).add(this.origin), kr));
      }
    }
    
    attractor_gen_poly(poly, n, kr) {
      let bounds = poly.get_bounds();
      for (let i = 0; i < n; i++) {
        let p = new dvector(dRandom.random(bounds[0].x, bounds[1].x), dRandom.random(bounds[0].y, bounds[1].y));
        //console.log(p)
        if (dGeom.point_in_poly(p, poly)) { 
          this.attractors.push(new dSpace_Growth_Attractor(p, kr));
        } else if (dGeom.point_in_poly(p.div(2), poly)) {
          this.attractors.push(new dSpace_Growth_Attractor(p.div(2), kr));
        } else {
          this.attractors.push(new dSpace_Growth_Attractor(p.normalize(), kr));
        }
      }
    }
  
    update(rc=false, r=1, poly=null) {
      let connections = {};
      let kill = [];
      for (let i = 0; i < this.attractors.length; i++) {
        let min = [0, 10000]
        for (let j = 0; j < this.nodes.length; j++) {
          let dist = this.nodes[j].pos.sub(this.attractors[i].pos).mag();
          if (dist < min[1] && (!rc || (rc && dist < r))) {
            min[0] = j;
            min[1] = dist;
          }
        }
        if (min[1] < this.attractors[i].kill_radius) { kill.push(this.attractors[i]); }
        if (!(min[0] in connections)) { connections[min[0]] = []; }
        connections[min[0]].push(this.attractors[i]);
  
      }
  
      for (const key_index in connections) {
        let value = connections[key_index];
        let key = this.nodes[key_index];
        //console.log(key, value);
        let av = new dvector(0, 0);
        for (let i = 0; i < value.length; i++) {
          av.add_equals(value[i].pos);
        }
        av.div_equals(value.length);
        let dir = av.sub(key.pos).normalize();
        let n = new dSpace_Growth_Node(key.pos.add(dir.mult(this.interval)), key);
        if (poly == null) { this.nodes.push(n); key.children.push(n); }
        else if (dGeom.point_in_poly(n.pos, poly)) { this.nodes.push(n); key.children.push(n); }
        
      }
  
      for (let i = 0; i < kill.length; i++) {
        let index = this.attractors.indexOf(kill[i]);
        this.attractors.splice(index, 1);
      }
  
      return connections;
    }
  }
}