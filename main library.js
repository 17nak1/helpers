var mathLib = {}
var erf = require('math-erf')


mathLib.pnorm = function (x, mu = 0, sd = 1, lower_tail = true, give_log = false) {
  if (sd < 0) {
    return NaN
  }
  let ans = 1 / 2 * (1 + erf((x - mu) / sd / Math.sqrt(2)))
  if (!lower_tail) {
    ans = 1 - ans
  }
  if (give_log) {
    ans = Math.log(ans)
  }
  return ans
}

mathLib.rnorm = function (mu = 0, sd = 1) {
  var val = Math.sqrt(-2 * Math.log(Math.random())) * Math.cos(2 * Math.PI * Math.random())
  return val * sd + mu
}

mathLib.dpois = function (x, lambda) {
  let ans, total = 0
  if (isNaN(x) || isNaN(lambda) || lambda < 0) {
    return NaN
  }
  if (!Number.isInteger(x)) {
    return 0
  }
  if (x < 0 || !isFinite(x)) {
    return 0
  }
  x = Math.round(x)
  ans = -lambda + x * Math.log(lambda)
  for (let i = 1; i <= x; i++) {
    total += Math.log(i)
  }
  let logAns = ans - total
  return Math.exp(logAns)
}

sum = function (array) {
  var sum = []  
  for(i = 0; i < array[0].length; i++){
    var s= 0
    for (j = 0; j < array.length; j++) {
       s += array[j][i] 
    }
    sum.push(s)
  }
  return sum
}

sp = function (scalar, array) {
  var sum = []
  for(i = 0; i < array.length; i++){
   sum.push(scalar * array[i]);
  }
  return sum
}
abs = function (array) {
  var sum = 0
  for(i = 0; i < array.length; i++){
   sum += Math.pow(Math.abs(array[i]), 2)
  }
  return Math.sqrt(sum)
}


mathLib.ode = function (method, func, N, params, t, h, pop, birthrate) {
  var tempArray
  var k1, k2, k3, k4, k5, k6, y, z, s
  var a, b, b2, c, d, out
  var tol = 1e-5
  switch (method) {
  case 'euler':
    tempArray = func(params, t, N, pop, birthrate)
    return [sum ([N,sp(h,tempArray)]), h]
    
  case 'rk4':
    c = [0, 1/3, 2/3, 1]
    a21 = 1/3 ,a31 = -1/3 ,a32 = 1, a41 = 1 ,a42 = -1 , a43 = 1
    b = [0, 1/8, 3/8, 3/8, 1/8]   
    k1 = func(params, t          , N, pop, birthrate)
    k2 = func(params, t + c[2] * h , sum([N , sp(h * a21 , k1)]), pop, birthrate)
    k3 = func(params, t + c[3] * h , sum([N , sp(h * a31 , k1), sp(h * a32, k2)]), pop, birthrate)
    k4 = func(params, t + c[4] * h , sum([N , sp(h * a41, k1),  sp(h *a42, k2), sp(h *a43, k3)]), pop, birthrate)
    return [sum ([N, sp (h *  b[1] , sum ([k2 , k3])) ,sp(h * b[2] ,sum ([k1 , k2]))]), h]
    
  case 'rkf45':
    c = [0,0, 1/4, 3/8, 12/13, 1, 1/2]
    a = [[0, 0, 0, 0, 0],
             [0,0, 0, 0, 0, 0],
             [0,1/4, 0, 0, 0, 0],
             [0,3/32, 9/32, 0, 0, 0],
             [0,1932/2197, -7200/2197, 7296/2197, 0, 0],
             [0,439/216, -8, 3680/513, -845/4104, 0],
             [0,-8/27, 2, -3544/2565, 1859/4104, -11/40]]      
    b = [0,25/216, 0, 1408/2565, 2197/4104, -1/5, 0]
    b2 = [0,16/135,   0,  6656/12825,   28561/56430,  -9/50,  2/55]
    k1 = func(params, t          , N, pop, birthrate)
    k2 = func(params, t + c[2] * h , sum([N , sp(h * a[2][1] , k1)]), pop, birthrate)
    k3 = func(params, t + c[3] * h , sum([N , sp(h * a[3][1] , k1), sp(h * a[3][2], k2)]), pop, birthrate)
    k4 = func(params, t + c[4] * h , sum([N , sp(h * a[4][1], k1),  sp(h *a[4][2], k2), sp(h *a[4][3], k3)]), pop, birthrate)
    k5 = func(params, t + c[5] * h , sum([N , sp(h * a[5][1], k1), sp(h *a[5][2], k2), sp(h *a[5][3], k3), sp(h *a[5][4], k4)]), pop, birthrate)
    k6 = func(params, t + c[6] * h , sum([N , sp(h * a[6][1], k1), sp(h *a[6][2], k2), sp(h *a[6][3], k3), sp(h *a[6][4], k4), sp(h *a[6][5], k5)]), pop, birthrate)
    y = sum ([N, sp (h *  b[1], k1), sp (h * b[2], k2) ,sp(h * b[3], k3), sp (h *  b[4], k4), sp (h * b[5], k5) ,sp(h * b[6], k6)])
    z = sum ([N, sp (h *  b2[1], k1), sp (h * b2[2], k2) ,sp(h * b2[3], k3), sp (h *  b2[4], k4), sp (h * b2[5], k5) ,sp(h * b2[6], k6)])
    s = abs(sum([z,sp(-1,y)])) ?  Math.sqrt(Math.sqrt((tol * h)/ (2 * abs(sum([z,sp(-1,y)]))))) : h;
    console.log((tol * h), (2 * abs(sum([z,sp(-1,y)]))))
    return[z, s] 
    
  case 'rk45dp7':
    c = [0, 0, 1/5, 3/10, 4/5, 8/9, 1, 1]
    a = [[0,0, 0, 0, 0, 0, 0],
            [0,0, 0, 0, 0, 0, 0],
            [0, 1/5, 0, 0, 0, 0, 0],
            [0, 3/40, 9/40, 0, 0, 0, 0],
            [0, 44/45, -56/15, 32/9, 0, 0, 0],
            [0, 19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0],
            [0, 9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0],
            [0, 35/384, 0, 500/1113, 125/192, -2187/6784, 11/84]]      
    b = [0, 5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40]
    b2 = [0, 35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0]
    d = [-12715105075/11282082432, 0, 87487479700/32700410799, -10690763975/1880347072, 701980252875/199316789632, -1453857185/822651844, 69997945/29380423]  
    k1 = func(params, t          , N, pop, birthrate)
    k2 = func(params, t + c[2] * h , sum([N , sp(h * a[2][1] , k1)]), pop, birthrate)
    k3 = func(params, t + c[3] * h , sum([N , sp(h * a[3][1] , k1), sp(h * a[3][2], k2)]), pop, birthrate)
    k4 = func(params, t + c[4] * h , sum([N , sp(h * a[4][1], k1),  sp(h *a[4][2], k2), sp(h *a[4][3], k3)]), pop, birthrate)
    k5 = func(params, t + c[5] * h , sum([N , sp(h * a[5][1], k1), sp(h *a[5][2], k2), sp(h *a[5][3], k3), sp(h *a[5][4], k4)]), pop, birthrate)
    k6 = func(params, t + c[6] * h , sum([N , sp(h * a[6][1], k1), sp(h *a[6][2], k2), sp(h *a[6][3], k3), sp(h *a[6][4], k4), sp(h *a[6][5], k5)]), pop, birthrate)
    y = sum ([N, sp (h *  b[1], k1), sp (h * b[2], k2) ,sp(h * b[3], k3), sp (h *  b[4], k4), sp (h * b[5], k5) ,sp(h * b[6], k6)])
    z = sum ([N, sp (h *  b2[1], k1), sp (h * b2[2], k2) ,sp(h * b2[3], k3), sp (h *  b2[4], k4), sp (h * b2[5], k5) ,sp(h * b2[6], k6)])
    s = Math.sqrt(Math.sqrt((tol * h)/ (2 * abs(sum([z,sp(-1,y)])))));
    // console.log(h, (2 * abs(sum([z,sp(-1,y)]))))
    return [z, s]
    
  }
  // return out  
}

mathLib.interpolator = function (points) {
  var first, n = points.length - 1,
    interpolated,
    leftExtrapolated,
    rightExtrapolated;

  if (points.length === 0) {
    return function () {
      return 0
    }
  }

  if (points.length === 1) {
    return function () {
      return points[0][1]
    }
  }

  points = points.sort(function (a, b) {
    return a[0] - b[0]
  })
  first = points[0]

  leftExtrapolated = function (x) {
    var a = points[0], b = points[1];
    return a[1] + (x - a[0]) * (b[1] - a[1]) / (b[0] - a[0])
  }

  interpolated = function (x, a, b) {
    return a[1] + (x - a[0]) * (b[1] - a[1]) / (b[0] - a[0])
  }

  rightExtrapolated = function (x) {
    var a = points[n - 1], b = points[n];
    return b[1] + (x - b[0]) * (b[1] - a[1]) / (b[0] - a[0])
  }

  return function (x) {
    var i
    if (x <= first[0]) {
      return leftExtrapolated(x)
    }
    for (i = 0; i < n; i += 1) {
      if (x > points[i][0] && x <= points[i + 1][0]) {
        return interpolated(x, points[i], points[i + 1])
      }
    }
    return rightExtrapolated(x);
  }
}

mathLib.sign = function (x, signal) {
  if (isNaN(x))
      return x
  return signal ? Math.abs(x) : -Math.abs(x);
}

mathLib.expRand = function (uniformRand) {
    var q = [
        0.6931471805599453,
        0.9333736875190459,
        0.9888777961838675,
        0.9984959252914960,
        0.9998292811061389,
        0.9999833164100727,
        0.9999985691438767,
        0.9999998906925558,
        0.9999999924734159,
        0.9999999995283275,
        0.9999999999728814,
        0.9999999999985598,
        0.9999999999999289,
        0.9999999999999968,
        0.9999999999999999,
        1.0000000000000000
    ];
    var a = 0.;
    var u = uniformRand();
    while (u <= 0. || u >= 1.)
        u = uniformRand();
    while (true) {
        u += u;
        if (u > 1.)
            break;
        a += q[0];
    }
    u -= 1.;
    if (u <= q[0])
      return a + u
    var i = 0;
    var ustar = uniformRand();
    var umin = ustar;
    do {
        ustar = uniformRand();
        if (umin > ustar)
            umin = ustar;
        i++;
    } while (u > q[i]);
    return a + umin * q[0];
}
mathLib.normalRand =function() {
    var u = 0, v = 0;
    while(u === 0) u = Math.random(); //Converting [0,1) to (0,1)
    while(v === 0) v = Math.random();
    return Math.sqrt( -2.0 * Math.log( u ) ) * Math.cos( 2.0 * Math.PI * v )
}

mathLib.numEulerSteps = function(t1, t2, dt) {
  var DOUBLE_EPS = 10e-8
  var tol = Math.sqrt(DOUBLE_EPS)
  var nstep
  if (t1 >= t2) {
    dt = 0
    nstep = 0
  } else if (t1 + dt >= t2) {
    dt = t2 - t1
    nstep = 1
  } else {
    nstep = Math.ceil((t2 - t1) / dt /(1 + tol))
    dt = (t2 - t1) / nstep
  }
  return nstep
}

mathLib.numMapSteps = function (t1, t2, dt) {
  var DOUBLE_EPS = 10e-8
  var tol = Math.sqrt(DOUBLE_EPS)
  var nstep
  // nstep will be the number of discrete-time steps to take in going from t1 to t2.
  nstep = Math.floor((t2 - t1) / dt /(1 + tol))
  return (nstep > 0) ? nstep : 0
}

mathLib.nosortResamp = function (nw, w, np, p, offset) {
  for (j = 1; j < nw; j++) {
   w[j] += w[j-1]
 }
  if (w[nw - 1] <= 0) {
    throw "in 'systematic_resampling': non-positive sum of weight"
  }
  var du = w[nw - 1] / np
  var u = -du * U.unif_rand()//Math.random()// U.unif_rand()

  for (j = 0, i = 0; j < np; j++) {
    u += du
    while ((u > w[i]) && (i < nw - 1)) i++;//looking for the low weight
    p[j] = i
  }
  if (offset){// add offset if needed
    for (j = 0; j < np; j++) p[j] += offset
  }
}


mathLib.reulermultinom = function (m = 1, size, rateAdd, dt, transAdd, rate, trans) {
  var p = 0
  var j, k
  if ((size < 0) || (dt < 0) || (Math.floor(size + 0.5) !== size)) {
    for (k = 0; k < m; k++) trans[k + transAdd] = NaN
    return 0
  }
  for (k = 0; k < m; k++) {
    if (rate[k + rateAdd] < 0.0) {
      for (j = 0; j < m; j++) trans[j + transAdd] = NaN
      return 0
    }
    p += rate[k + rateAdd]// total event rate
  }
  if (p > 0) {
    size = rbinom.rbinomOne(size, 1 - Math.exp(-p * dt))// total number of events
    if (!(isFinite(size)))
      throw 'result of binomial draw is not finite.'
    m -= 1
    for (k = 0; k < m; k++) {
      if (rate[k + rateAdd] > p) p = rate[k + rateAdd]
      trans[k + transAdd] = ((size > 0) && (p > 0)) ? rbinom.rbinomOne(size, rate[k + rateAdd] / p) : 0
      if (!(isFinite(size) && isFinite(p) && isFinite(rate[k + rateAdd]) && isFinite(trans[k + transAdd]))) {
        throw 'result of binomial draw is not finite.'
      }
      size -= trans[k + transAdd]
      p -= rate[k + rateAdd]
    }
    trans[m + transAdd] = size
  } else {
    for (k = 0; k < m; k++) trans[k + transAdd] = 0
  }
}

mathLib.logMeanExp = function (x) {
  var mx = Math.max(...x)
  var s = x.map((x,i) => Math.exp(x - mx))
  var q = s.reduce((a, b) => a + b, 0)
  return mx + Math.log(q / x.length)
}
module.exports = mathLib



