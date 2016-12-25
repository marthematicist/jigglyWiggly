// waveMesh
// marthematicist - 2016
var vers = '0.21';
console.log( 'jigglyWiggly: version ' + vers );


// GLOBAL VARIABLES /////////////////////////////////////////
function setupGlobalVariables() {
  // DISPLAY WINDOW VARIABLES
  {
    // window resolution (pixels)
    xRes = windowWidth;
    yRes = windowHeight;
    minRes = min( xRes , yRes );
    maxRes = max( xRes , yRes );
    winArea = xRes*yRes;
  }
  
  // DRAW VARIABLES
  {
    // bachground transparency
    bgAlpha = 255;
    // background color
    bgColor = color( 255 , 255 , 255 , bgAlpha );
    // node transparency
    nodeAlpha = 255;
    // node color
    nodeColor = color( 255 , 255 , 255 , nodeAlpha );
    // node diameter
    nodeDiam = 0.01*minRes;
    // mesh transparency
    meshAlpha = 255;
    // mesh color
    meshColor = color( 0 , 0 , 0 , meshAlpha );
    // mesh thickness
    meshWeight = 0.003*minRes;
    // fill transparency
    fillAlpha =255;
    // fill base colors
    minc = 150;
    maxc = 255;
    minColorDev = 20;
    fillBaseColor1 = randomColor( minc , maxc , minColorDev , fillAlpha );
    //fillBaseColor1 = color( random(minc,maxc) , random(minc,maxc) , random(minc,maxc) , fillAlpha );
    fillBaseColor2 = randomColor( minc , maxc , minColorDev , fillAlpha );
    // fill color random deviation
    minFillLerp = 0.1;
    maxFillLerp = 0.2;
  }
  
  // SIMULATION VARIABLES
  {
    // desired number of nodes
    desiredNodes = 200;
    // total mass
    totalMass = desiredNodes*2;
    // simulation area
		simArea = desiredNodes;
		// linear conversion factor: sim to window
		sim2WinFactor = sqrt( winArea / simArea );
		// linear conversion factor: window to sim
		win2SimFactor = sqrt( simArea / winArea );
		// dimensions of the simulation
		xExt = xRes * win2SimFactor;
		yExt = yRes * win2SimFactor;
		// the minimum dimension of the simulation
		minExt = min( xExt , yExt );
		// bounds of simulation
		xMin = -0.5*xExt;
		yMin = -0.5*yExt;
		xMax = 0.5*xExt;
		yMax = 0.5*yExt;
		// simulation center (vector)
		simCenter = createVector( 0.5*( xMin + xMax ) , 0.5*( yMin + yMax ) );
		// mesh dimensions
		N = round( yExt );
		M = round( xExt );
		numNodes = M*N;
		// node separation
		xSep = xExt / ( N - 1 );
		ySep = yExt / ( M - 1 );
		// node mass
		massPerNode = totalMass / numNodes;
		// perturber mass
    perturberMass = 3*massPerNode;
  }
  
  // PHYSICS VARIABLES
  {
    // perturber repulse?
    perturberRepulse = false;
    // max perturber force for attract
    maxPerturberForceAttract = 0.25;
    // factor relating simulation time to frame time
    // frameTime*dtt = simTime
    dtt = 0.20 / 20;
    // simulation time per frame
		dt = 0.20;
		// gravity constant
		gravityConstant = 1;
		// gravity smoothing factor (0.4 optimum)
		epsilon = 2;
		// spring constant
		springConstant = 1;
		// friction coefficient
		frictionCoefficient = 0.05;
		// whether or not perturber exists
		perturberExists = false;
  }
  
  // TIMING VARIABLES
  {
    // time that setup was completed
    startTimer = 0;
		// amount of time to wait after setup before beginning simulation (title screen time)
		startWaitTime = 4000;
		// used to clear the screen after title
		clearFirstTime = true;
		// time since last color change
		colorTimer = 0;
		// time between color changes
		colorWaitTime = 10000;
		// time of last cursor move
		moveTimer = millis();
		// delay between last move and deactivation of perturber
		moveWaitTime = 600;
  }
  
  //RECORD-KEEOING VARIABLES
  {
    // previous values of mouseX, mouseY
    lastMouseX = mouseX;
    lastMouseY = mouseY;
    // milliseconds to draw a frame
    frameTimer = millis();
    // rolling averave of frameTimer
    avgFrameTime = 20;
    // maximum force magnitude from perturber
    maxfMag = 0;
  }
  // GLOBAL OBJECTS
  {
    mesh = false;
  }
  
  // DISPLAY SETUP VALUES TO CONSOLE
  {
    console.log( 'M=' + M + '   N=' + N + '   num=' + numNodes );
  }
}

// UTILITY FUNCTIONS //////////////////////////////////////////////////////
// functions to convert linear distance and vectors from sim to window
function sim2Win( a ) {
	return a*sim2WinFactor;
}
function sim2WinVect( a ) {
	return createVector( (a.x - xMin)*sim2WinFactor ,
						 (a.y - yMin)*sim2WinFactor );
}
function win2SimVect( a ) {
	return createVector( (a.x - 0.5*xRes)*win2SimFactor ,
						 (a.y - 0.5*yRes)*win2SimFactor);
}
// function to choose a random color
function randomColor( minCV , maxCV , minDev , alpha ) {
  var colorAcceptable = false;
  var r = 0;
  var g = 0;
  var b = 0;
  while( !colorAcceptable ) {
    r = random( minCV , maxCV );
    g = random( minCV , maxCV );
    b = random( minCV , maxCV );
    var a = (r + g + b)/3;
    var dev = sqrt( ( (a-r)*(a-r) + (a-g)*(a-g) + (a-b)*(a-b) ) / 3 );
    if( dev > minDev ) { colorAcceptable = true; }
  }
  return color( r , g , b , alpha );
}

// CLASS: Node /////////////////////////////////////////////////////////////
var Node = function( m , n , x , mass ) {
  // CONSTRUCTOR INPUTS:
  //  m : vertical position in mesh
  //  n : horizontal position in mesh
  //  x : physical position (2D p5.Vector)
  //  mass : mass of node
  
  // CLASS VARIABLES:
  // vertical (m) and horizontal positions in mesh
  this.m = m;
  this.n = n;
  // position in simulation space
  this.x = createVector( x.x , x.y );
  // velocity vector
  this.v = createVector( 0 , 0 );
  // acceleration vector
  this.a = createVector( 0 , 0 );
  // mass
  this.mass = mass;
  // array of Boolean: does it have neighbors? (set in Mesh method)
  // [0] right ; [1] top ; [2] left ; [3] bottom
  this.hasNeighbors = new Array( 4 );
  // array of Nodes: the neighbors (set in Mesh method)
  this.neighbors = new Array( 4 );
  // is this node's position fixed? (set border nodes in Mesh method)
  this.fixed = false;
  // node color
  this.color = color( random(minc,maxc) , random(minc,maxc) , random(minc,maxc) , fillAlpha );
  // lerp amount
  this.lerpAmt = random( minFillLerp , maxFillLerp );
};

// CLASS: Mesh /////////////////////////////////////////////////////////////
var Mesh = function( MM , NN , xmin , xmax , ymin , ymax , tMass ) {
  // CONSTRUCTOR INPUTS:
  //  MM : number of rows
  //  NN : number of columns
  //  xmin , xmax , ymin , ymax : mesh edges
  //  tMass : total mass
  
  // CLASS VARIABLES:
  // number of rows
  this.M = MM;
  // number of columns
  this.N = NN;
  // number of nodes
  this.num = this.N * this.M;
  // edges and extent
  this.xMin = xmin;
  this.xMax = xmax;
  this.yMin = ymin;
  this.yMax = ymax;
  this.xExt = this.xMax - this.xMin;
  this.yExt = this.yMax - this.yMin;
  this.minExt = min( this.xExt , this.yExt );
  this.maxExt = max( this.xExt , this.yExt );
  // total mass
  this.totalMass = tMass;
  // array of Node objects
  this.nodes = new Array( this.num );
  var dx = this.xExt / (this.M - 1);
  var dy = this.yExt / (this.N - 1);
  var mass = this.totalMass / this.num;
  for( var m = 0 ; m < this.M ; m++ ) {
    for( var n = 0 ; n < this.N ; n++ ) {
      var ind = m*this.N + n;
      var x = createVector( this.xMin + m*dx , this.yMin + n*dy );
      this.nodes[ind] = new Node( m , n , x , mass );
    }
  }
  // perturber position
  this.px = createVector( 0 , 0 );
  // perturber mass
  this.pmass = perturberMass;

  
  // METHODS:
  // Mesh method: node
  // returns Node with indices m,n
  this.node = function( m , n ) {
    return( this.nodes[ m*this.N + n ] );
  };
  
  // Mesh method: setNeighbors
  // sets neighbor info for all nodes, and sets border nodes as fixed
  this.setNeighbors = function() {
    // for each node...
    for( var i = 0 ; i < this.num ; i++ ) {
      var m = this.nodes[i].m;
      var n = this.nodes[i].n;
      // check/set right neighbor
      if( m === this.M-1 ) {
        // if node is on right border, set right neighbor to false,
        // and fix the node's position
        this.nodes[i].hasNeighbors[0] = false;
        this.nodes[i].neighbors[0] = false;
        this.nodes[i].fixed = true;
      } else {
        // otherwise the node has a neighbor on the right and is not fixed
        this.nodes[i].hasNeighbors[0] = true;
        this.nodes[i].neighbors[0] = this.node( m+1 , n );
      }
      // check/set top neighbor
      if( n === 0 ) {
        // if node is on top border, set top neighbor to false,
        // and fix the node's position
        this.nodes[i].hasNeighbors[1] = false;
        this.nodes[i].neighbors[1] = false;
        this.nodes[i].fixed = true;
      } else {
        // otherwise the node has a neighbor on the top and is not fixed
        this.nodes[i].hasNeighbors[1] = true;
        this.nodes[i].neighbors[1] = this.node( m , n-1 );
      }
      // check/set left neighbor
      if( m === 0 ) {
        // if node is on left border, set left neighbor to false,
        // and fix the node's position
        this.nodes[i].hasNeighbors[2] = false;
        this.nodes[i].neighbors[2] = false;
        this.nodes[i].fixed = true;
      } else {
        // otherwise the node has a neighbor on the left and is not fixed
        this.nodes[i].hasNeighbors[2] = true;
        this.nodes[i].neighbors[2] = this.node( m-1 , n );
      }
      // check/set bottom neighbor
      if( n === this.N-1 ) {
        // if node is on bottom border, set bottom neighbor to false,
        // and fix the node's position
        this.nodes[i].hasNeighbors[3] = false;
        this.nodes[i].neighbors[3] = false;
        this.nodes[i].fixed = true;
      } else {
        // otherwise the node has a neighbor on the bottom and is not fixed
        this.nodes[i].hasNeighbors[3] = true;
        this.nodes[i].neighbors[3] = this.node( m , n+1 );
      }
    }
  };
  
  // Mesh method: drawNodes
  // draws nodes to canvas
  this.drawNodes = function() {
    noStroke();
    fill( nodeColor );
    for( var i = 0 ; i < this.num ; i++ ) {
      var x = sim2WinVect( this.nodes[i].x );
      ellipse( x.x , x.y , nodeDiam , nodeDiam );
    }
  };
  
  
  // Mesh method: drawMesh
  // draws lines between neighboring nodes
  this.drawMesh = function() {
    stroke( meshColor );
    strokeWeight( meshWeight );
    noFill();
    for( var m = 0 ; m < this.M ; m++ ) {
      beginShape();
      for( var n = 0 ; n < this.N ; n++ ) {
        var v = sim2WinVect( this.node(m,n).x );
        vertex( v.x , v.y );
      }
      endShape();
    }
    for( var n = 0 ; n < this.N ; n++ ) {
      beginShape();
      for( var m = 0 ; m < this.M ; m++ ) {
        var v = sim2WinVect( this.node(m,n).x );
        vertex( v.x , v.y );
      }
      endShape();
    }
  };
  
  // Mesh method: drawMeshFill
  // fills in the spaces in the mesh
  this.drawMeshFill = function( baseColor ) {
    noStroke();
    for( m = 0 ; m < this.M - 1 ; m++ ) {
      for( n = 0 ; n < this.N - 1 ; n++ ) {
        var fc =
        fill( lerpColor( baseColor , this.node(m,n).color , this.node(m,n).lerpAmt ) );
        var v0 = sim2WinVect( this.node(m,n).x );
        var v1 = sim2WinVect( this.node(m+1,n).x );
        var v2 = sim2WinVect( this.node(m+1,n+1).x );
        var v3 = sim2WinVect( this.node(m,n+1).x );
        beginShape();
        vertex( v0.x , v0.y );
        vertex( v1.x , v1.y );
        vertex( v2.x , v2.y );
        vertex( v3.x , v3.y );
        endShape( CLOSE );
      }
    }
  };
  
  // Mesh method: zeroAccelerations
  this.zeroAccelerations = function() {
    for( var i = 0 ; i < this.num ; i++ ) {
      this.nodes[i].a = createVector( 0 , 0 );
    }
  };
  
  // Mesh method: applyFrictionForces
  this.applyFrictionForces = function() {
    for( var i = 0 ; i< this.num ; i++ ) {
      if( !this.nodes[i].fixed ) {
        var dA = createVector( this.nodes[i].v.x , this.nodes[i].v.y );
        dA.mult( -frictionCoefficient );
        this.nodes[i].a.add( dA );
      }
    }
  }
  
  // Mesh method: applySpringForces
  this.applySpringForces = function() {
    for( var i = 0 ; i< this.num ; i++ ) {
      var x1 = this.nodes[i].x;
      var m1 = this.nodes[i].mass;
      // if node has right neighbor...
      if( this.nodes[i].hasNeighbors[0] ) {
        // get the neighbor position and mass
        var x2 = this.nodes[i].neighbors[0].x;
        var m2 = this.nodes[i].neighbors[0].mass;
        // find the direction and magnitude from 1 to 2
        var fd = p5.Vector.sub( x2 , x1 );
        var fm = fd.mag() * springConstant;
        fd.normalize();
        // compute force vector (1 to 2)
        fd.mult( fm );
        if( !this.nodes[i].fixed ) {
          this.nodes[i].a.add( p5.Vector.div( fd , m1 ) );
        }
        if( !this.nodes[i].neighbors[0].fixed ) {
          this.nodes[i].neighbors[0].a.sub( p5.Vector.div( fd , m2 ) );
        }
      }
      // if node has top neighbor...
      if( this.nodes[i].hasNeighbors[1] ) {
        // get the neighbor position and mass
        var x2 = this.nodes[i].neighbors[1].x;
        var m2 = this.nodes[i].neighbors[1].mass;
        // find the direction and magnitude from 1 to 2
        var fd = p5.Vector.sub( x2 , x1 );
        var fm = fd.mag() * springConstant;
        fd.normalize();
        // compute force vector (1 to 2)
        fd.mult( fm );
        if( !this.nodes[i].fixed ) {
          this.nodes[i].a.add( p5.Vector.div( fd , m1 ) );
        }
        if( !this.nodes[i].neighbors[1].fixed ) {
          this.nodes[i].neighbors[1].a.sub( p5.Vector.div( fd , m2 ) );
        }
      }
    }
  }
  
  // Mesh method: applyPerturberForces
  this.applyPerturberForces = function() {
    for( var i = 0 ; i< this.num ; i++ ) {
      if( !this.nodes[i].fixed ) {
        var d = p5.Vector.dist( this.px , this.nodes[i].x );
        var acc = p5.Vector.sub( this.px , this.nodes[i].x );
        acc.normalize();
        var fMag = this.pmass * gravityConstant / pow( d * d + epsilon * epsilon , 1.5 );
        acc.mult( fMag );
        this.nodes[i].a.sub( acc );
      }
    }
  };
  
  // Mesh method: evolveHalfStep
  // evolves the mesh physics 1/2 step (leapfrog integrator)
  this.evolveHalfStep = function() {
    // zero, then compute new accelerations
    this.zeroAccelerations();
    if( perturberExists ) {
      this.applyPerturberForces();
    }
    this.applySpringForces();
    this.applyFrictionForces();
    // update velocities (1/2 step)
    for( var i = 0 ; i < this.num ; i++ ) {
      this.nodes[i].v.add( p5.Vector.mult( this.nodes[i].a , dt*0.5 ) );
    }
  };
  
  // Mesh method: evolveFullStep
  // evolves the mesh physics (steps) full steps (leapfrog integrator)
  this.evolveFullStep = function( steps ) {
    for( var s = 0 ; s < steps ; s++ ) {
      // update positions
      for( var i = 0 ; i < this.num ; i++ ) {
        this.nodes[i].x.add( p5.Vector.mult( this.nodes[i].v , dt ) );
      }
      // zero, then compute new accelerations
      this.zeroAccelerations();
      if( perturberExists ) {
        this.applyPerturberForces();
      }
      this.applySpringForces();
      this.applyFrictionForces();
      // update velocities (full step)
      for( var i = 0 ; i < this.num ; i++ ) {
        this.nodes[i].v.add( p5.Vector.mult( this.nodes[i].a , dt ) );
      }
    }
  };
};
// function to create and initialize a new mesh
function createMesh( MM , NN , xmin , xmax , ymin , ymax , tMass ) {
  var output = new Mesh( MM , NN , xmin , xmax , ymin , ymax , tMass );
  output.setNeighbors();
  return output;
}

// p5 SETUP FUNCTION ///////////////////////////////////////////////////////
function setup() {
  // set up global variables
  setupGlobalVariables();
  // define canvas
  createCanvas( xRes , yRes );
  // initialize Mesh
  mesh = createMesh( M , N , xMin , xMax , yMin , yMax , totalMass );
  // set base fill colors
  fillBaseColor1 = color( random(minc,maxc) , random(minc,maxc) , random(minc,maxc) , fillAlpha );
  fillBaseColor2 = color( random(minc,maxc) , random(minc,maxc) , random(minc,maxc) , fillAlpha );
  // evolve the mesh 1/2 step
  mesh.evolveHalfStep();
  
  startTimer = millis();
}

// p5 DRAW FUNCTION /////////////////////////////////////////////////////////
function draw() {
  
  // if mouse has moved...
  if( lastMouseX != mouseX || lastMouseY != mouseY ) {
    // reset move timer
    moveTimer = millis();
    // turn on the pertueber
    perturberExists = true;
    mesh.pmass = perturberMass;
  } else {
    // if mouse hasn't moved, decrease perturber mass
    var t = millis() - moveTimer;
    if( t > moveWaitTime ) {
      // if move wait time exceeded, turn off the perturber
      perturberExists = false;
    } else {
      mesh.pmass = ( 1 - t/moveWaitTime )*perturberMass;
    }
  }
  // set lastMouseX and lastMouseY
  lastMouseX = mouseX;
  lastMouseY = mouseY;
  
  // if color change time exceeded, change colors
  if( (millis() - colorTimer) > colorWaitTime ) {
    // set base fill 1 to previous base fill 2
    fillBaseColor1 = color( red(fillBaseColor2) , green(fillBaseColor2) ,
                            blue(fillBaseColor2) , fillAlpha );
    // set base fill 2 to a random color
    fillBaseColor2 = randomColor( minc , maxc , minColorDev , fillAlpha );
    // reset color timer
    colorTimer = millis();
  }
  
  // set base color
  var lerpAmt = (millis() - colorTimer) / colorWaitTime;
  var fillBaseColor = lerpColor( fillBaseColor1 , fillBaseColor2 , lerpAmt );
  
  // set perturber to mouse position
  var p = win2SimVect( createVector( mouseX , mouseY ) );
  mesh.px = createVector( p.x , p.y );
  
  // evlolve physics full step
  mesh.evolveFullStep( 1 );
  
  // draw background
  background( bgColor );
  // draw nodes and mesh and fill
  mesh.drawMeshFill( fillBaseColor );
  mesh.drawMesh();
  //mesh.drawNodes();
  
  // if still in splash display, draw splash over screen
	if( millis() - startTimer < startWaitTime ) {
		// display beginning text
  	background( 0 , 0 , 0 );
  	textAlign( CENTER );
  	textSize( minRes*0.08 );
  	fill(255);
  	text("JIGGLY-WIGGLY\n-marthematicist-" , 0.5*xRes , 0.5*yRes - 80 );
  	textSize( minRes*0.05 );
  	//text( "A physics simulation with springs" , 0.5*xRes , 0.5*yRes + 70 );
  	//textSize( minRes*0.04 );
  	text( "version " + vers , 0.5*xRes , yRes - 100 );
	} else if( clearFirstTime ) {
	  clearFirstTime = false;
	  dt = dtt * avgFrameTime;
    if( dt > 0.75 ) { dt = 0.75; }
    console.log( 'average frame time:' + avgFrameTime +'   dt=' + dt );
	}
  
  // update avgFrameTime
  avgFrameTime = 0.9*avgFrameTime + 0.1*( millis() - frameTimer );
  
  // reset frameTimer
  frameTimer = millis();
  
  // log out data periodically
  if( frameCount % 100 === 0 ) {
    // reset dt
    dt = dtt * avgFrameTime;
    if( dt > 0.75 ) { dt = 0.75; }
    console.log( 'average frame time:' + avgFrameTime +'   dt=' + dt );
    
  }
  //console.log( mesh.nodes[0].v );
}


