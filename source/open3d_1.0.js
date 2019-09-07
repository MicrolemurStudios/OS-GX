//This is a Processing.JS File. It needs the Processing.JS Libary to run!

/**
Open3D is a open-source libary for 3D Games. It is still very early in development and in alpha stages.

Version History:
 
Open3D 1.0:

Basic Cubes/Coins Ready for action now.

*/
var open3d_version="1.0";
/* RMatrix2D is copied from https://www.khanacademy.org/cs/mr/3793059675 */
var RMatrix2D = (function() {
    var aTan2 = atan2;  /* copy to a smaller, nearer, global name space. */
    var abs = Math.abs;
    var sqrt = Math.sqrt;
    var orthogEpsilon = 1/1024/16;
    var epsilon = 1/1024/1024/1024;
    var freeMatrices = [];  /* To mitigate KA storage leak... */
    
    /* Constructor function. Objects inherrit everything from PMatrix2D */
    var RMatrix2D = function() {
        PMatrix2D.apply(this, arguments);
        this.mops = {};
    };
    
    /*
     * Inherited PMatrix2D properties include:
     *  elements: a six element array of numbers
     *  set: function takes another PMatrix2D, an array, or six numbers
     *  get: function replaced by RMatrix2D.prototype.get
     *  reset: function sets this to the identity matrix
     *  translate: function
     *  determinant: function
     *  invert: function this = inverse of this, returns a success boolean
     *  scale: function
     *  apply: function this = this * argument
     *  preApply: function this = argument * this
     *  rotate: function
     *  print: function like printMatrix
     * Others that should most likely be avoided are:
     *  array: function returns a new array of this.elements
     *  mult: function
     *  multX: function
     *  multY: function
     *  skewX: function
     *  skewY: function
     *  transpose: function
     *  rotateZ: invokes rotate
     *  invRotateZ: function() undo rotateZ
     *  invTranslate: function undo translate
     *  invScale: function undo scale
     */
    RMatrix2D.prototype = Object.create(PMatrix2D.prototype);
    RMatrix2D.prototype.resetMatrix = PMatrix2D.prototype.reset;  /* alias */
    RMatrix2D.prototype.printMatrix = PMatrix2D.prototype.print;  /* alias */
    
    /* Use this instead of "new RMatrix2D" to reuse recycled matrices. */
    RMatrix2D.new = function() {
        var m = freeMatrices.pop() || new RMatrix2D();
        m.set.apply(m, arguments);
        return m;
    };
    
    /* Free/recycle this matrix when it is no longer needed. */
    RMatrix2D.prototype.free = function() {
        var stack = this.stack;
        while (stack && stack.length > 0) {
            stack.pop().free();
        }
        if (this.inverted) {
            this.inverted.free();
            this.inverted = null;
        }
        if (! this.breakMirror) {
            this.reset();
            freeMatrices.push(this);
        }
        return null;
    };
    
    /* Return a clone of this. */
    RMatrix2D.prototype.get = function() {
        return RMatrix2D.new(this);
    };
    
    /* Provide the same functionality as pushMatrix() as a method. */
    RMatrix2D.prototype.push = RMatrix2D.prototype.pushMatrix = function() {
        this.stack = this.stack || [];
        this.stack.push(RMatrix2D.new(this));
    };
    
    /* Provide the same functionality as popMatrix() as a method. */
    RMatrix2D.prototype.pop = RMatrix2D.prototype.popMatrix = function() {
        var stack = this.stack;
        if (stack && stack.length > 0) {
            var top = stack.pop();
            this.set(top);
            top.free();
        }
    };
    
    /*
     * Compute this matrix A times the PVector v.
     * The result PVector is set. results can be v.
     */
    RMatrix2D.prototype.transform = function(v, results) {
        var A = this.elements;
        results.set(A[0]*v.x + A[1]*v.y + A[2], A[3]*v.x + A[4]*v.y + A[5], 0);
        return results;
    };
    
    /*
     * Return the inverse of this matrix.  Caller provides
     * a target results matrix or is provided a STATIC
     * property of this matrix.
     */
    RMatrix2D.prototype.inverse = function(results) {
        results = results || (this.inverted = this.inverted || RMatrix2D.new());
        results.set(this);
        if (! results.invert()) {
            throw "No matrix inverse possible!";
        }
        return results;
    };
    
    /* Apply shear offsets to this matrix. */
    var shearMatrix = new RMatrix2D();
    RMatrix2D.prototype.shear = function(shx, shy) {
        shearMatrix.set(1, shx, 0, shy, 1, 0);
        this.apply(shearMatrix);
    };
    
    /* Shear the Processing.js rendering matrix. */
    RMatrix2D.shear = function(shx, shy) {
        shearMatrix.set(1, shx, 0, shy, 1, 0);
        shearMatrix.applyMopsToPJS();
    };
    
    /* 
     * Compute this 2×2 matrix's eigenvalues and UNIT eigenvectors. See
     * http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/index.html
     * for formulas.
     */
    RMatrix2D.prototype.eigen = function() {
        var ev1 = this.eigenVector1 = this.eigenVector1 || new PVector();
        var ev2 = this.eigenVector2 = this.eigenVector2 || new PVector();
        var E = this.elements;
        var trace = E[0] + E[4];
        var discrim = trace * trace / 4 - this.determinant();
        discrim = (discrim < epsilon) ? 0 : sqrt(discrim);
        var lamda1 = this.lamda1 = trace/2 + discrim;
        var lamda2 = this.lamda2 = trace/2 - discrim;
        if (abs(E[3]) > epsilon) {
            ev1.set(lamda1 - E[4], E[3], 0);
            ev2.set(lamda2 - E[4], E[3], 0);
        } else if (abs(E[1]) > epsilon) {
            ev1.set(E[1], lamda1 - E[0], 0);
            ev2.set(E[1], lamda2 - E[0], 0);
        } else if (abs(E[0] - lamda1) < epsilon) {
            ev1.set(1, 0, 0);
            ev2.set(0, 1, 0);
            return;
        } else {
            ev2.set(1, 0, 0);
            ev1.set(0, 1, 0);
            return;
        }
        ev1.normalize();
        ev2.normalize();
    };
    
    var myV = new RMatrix2D();
    var myP = new RMatrix2D();
    /*
     * Find the Singular Value Decomposition of this. Parameters
     * are all required RMatrix2D output targets. So this = U×S×VT.
     * Returns true on success. See section 8.1 of
     * https://datajobs.com/data-science-repo/SVD-Tutorial-[Kirk-Baker].pdf
     */
    RMatrix2D.prototype.svd = function(U, S, VT) {
        var E = this.elements, a = E[0], b = E[1], c = E[3], d = E[4];
        /* Degenerate matrices (collinear vectors) are a big fail. */
        if (abs(a*d - b*c) < epsilon) {
            return false;
        }
        
        /* Special (simple) case for orthogonal vectors. */
        if (abs(a*b + c*d) < orthogEpsilon) {
            var magX = sqrt(a*a + c*c);
            var magY = sqrt(b*b + d*d);
            U.set(  a/magX, b/magY, 0,
                    c/magX, d/magY, 0 );
            S.set(  magX, 0, 0,
                    0, magY, 0);
            VT.resetMatrix();  /* VT = Identity matrix */
            return true;
        }
        
        /* U is derived from the eigenvector of EE' */
        myP.set(a, c, 0, b, d, 0);  /* sets myP = E' */
        myP.preApply(this);
        myP.eigen();
        var u1 = myP.eigenVector1, u2 = myP.eigenVector2;
        U.set(u1.x, u2.x, 0, u1.y, u2.y, 0);
        
        /* S diag is 1/sqrt of the eigenvalues */
        S.set(sqrt(myP.lamda1), 0, 0, 0, sqrt(myP.lamda2), 0);
        
        /* 
         * V = [v1, v2] where the column vectors
         * v1 = (1/sqrt(lamda1))·E'·u1  and
         * v2 = (1/sqrt(lamda2))·E'·u2.
         */
        var V = myV;
        myP.set(a, c, 0, b, d, 0);  /* myP = E' */
        myP.transform(u1, u1);
        u1.div(S.elements[0]);
        myP.transform(u2, u2);
        u2.div(S.elements[4]);
        V.set(u1.x, u2.x, 0, u1.y, u2.y, 0);
        
        /*
         * Multiply U and V by minusOne if U does not appear to be a
         * valid rotation matrix.  minusOne × minusOne === Identity, so
         * U×S×V remains the same. minusOne = [ -1, 0, 0, 0, 1, 0 ].  This
         * is the same as negating each X column vector of the two matrices.
         */
        U = U.elements;
        V = V.elements;
        if (U[0] * U[4] <= 0) {
            U[0] *= -1;
            U[3] *= -1;
            V[0] *= -1;
            V[3] *= -1;
        }
        
        VT.set(V[0], V[3], 0, V[1], V[4], 0);  /* VT = V' */
        return true;
    };
    
    var myU = new RMatrix2D();
    var myS = new RMatrix2D();
    var myVT = new RMatrix2D();
    var rotReflect =  new RMatrix2D();
    var tuY = new PVector();
    /*
     * This method provides the arguments to translate,
     * rotate, and scale that can reconstruct a rendering
     * matrix. Results parameter is returned. If no results
     * parameter is given then a STATIC object is returned.
     * Will return null on failure.  See applyMopsToPJS()
     * for HOW to apply the transformations.
     */
    RMatrix2D.prototype.reconstitute = function(results) {
        var matrix = this;
        var mops = results || matrix.mops;
        var E = matrix.elements, a = E[0], b = E[1], c = E[3], d = E[4];
        mops.tx = E[2];
        mops.ty = E[5];
        if (a*d - b*c < 0) {
            /* 
             * A negative Z component of tuX ⨯ tuY means that the
             * Y axis is a negative angle from X axis in a universe 
             * where the Y axis must be a positive angle from X axis.
             * So, find the SVD of a similar "positive" matrix whose
             * X column vector is collinear with the X axis, and
             * whose Y column vector is reflected by the X axis.
             */
            tuY.set(b, d, 0);
            mops.rot0 = aTan2(c, a);  /* rotation FROM X axis */
            tuY.rotate(-mops.rot0);  /* normalize tuX & tuV back TO X axis. */
            rotReflect.set(sqrt(a*a + c*c), tuY.x, 0, 0, -tuY.y, 0);
            matrix = rotReflect;
            mops.initSy = -1;  /* reconcile reflection, -tuY.y */
        } else {
            mops.rot0 = 0;  /* no rotation */
            mops.initSy = 1;  /* no reflection */
        }
        
        /* Now, factor matrix into three matrices... */
        if (! matrix.svd(myU, myS, myVT)) {
            return null;
        }
        
        /*
         * Since matrix = U×S×VT, compute arguments to rotate and
         * scale that can reconstruct matrix. Do the first rotate...
         */
        mops.rot1 = aTan2(myU.elements[3], myU.elements[0]);
        /* Then do a scale. */
        mops.sx = myS.elements[0];
        mops.sy = myS.elements[4];
        /* Finally do another rotate. */
        mops.rot2 = aTan2(myVT.elements[3], myVT.elements[0]);
        return mops;
    };
    
    /*
     * Static function applies the results of the reconstitute()
     * method to THE Processing.js rendering matrix. The only
     * way to do that is via a series of tranformation calls.
     * Returns true if mops.
     */
    var applyMopsToPJS = RMatrix2D.applyMopsToPJS = function(mops) {
        if (mops) {
            void((mops.tx || mops.ty) && translate(mops.tx, mops.ty));
            void(mops.rot0 && rotate(mops.rot0));
            void((mops.initSy !== 1) && scale(1, mops.initSy));
            void(mops.rot1 && rotate(mops.rot1));
            void((mops.sx !== 1 || mops.sy !== 1) && scale(mops.sx, mops.sy));
            void(mops.rot2 && rotate(mops.rot2));
            return true;
        }
    };
    
    /*
     * Apply this RMatrix2D to THE Processing.js rendering matrix.
     * Returns true iff successful.
     */
    RMatrix2D.prototype.applyMopsToPJS = function() {
        return applyMopsToPJS(this.reconstitute());
    };
    
    /*
     * (The Holy Grail.)
     * Static function applies the RMatrix2D rm to THE 
     * Processing.js rendering matrix. If you want to SET
     * the Processing.js matrix, then invoke resetMatrix()
     * prior to invoking this function.  Returns true iff
     * successful.
     */
    RMatrix2D.applyMatrix = function(rm) {
        return rm.applyMopsToPJS();
    };
    
    /*
     * Static function that converts the string
     * produced by Processing.js printMatrix to a
     * RMatrix2D object.  KA printMatrix() outputs
     * six numbers on two lines that look like 
     *      a b c
     *      d e f
     * You should pass ONE string that looks like
     *      a b c d e f
     * where the numbers are separated by blanks.
     * Optional parameter results is where the result
     * are placed and returned.  It must be a RMatrix2.
     * If none is provided then a new one is created.
     */
    RMatrix2D.string2RMatrix = function(s, results) {
        try {  /* more KA silliness */
            results = results || RMatrix2D.new();
            var token = s.split(" ");
            var num = results.elements;
            for (var n = 0, i = 0; (n < 6) && (i < token.length); i++) {
                var f = parseFloat(token[i]);
                if (f > f - 1) {
                    /* f is really a number */
                    num[n++] = f;
                }
            }
            if (n !== 6) {
                println("Expected six, but got " + n + " numbers!");
            }
            return results;
        } catch(err) {}
    };
    
    /* Names of extremely interesting Processing.js functions. */
    RMatrix2D.transforms = [ "translate", "scale", "rotate",
        "pushMatrix", "popMatrix", "resetMatrix" ];
        
    /*
     * This static function returns a RMatrix2D that 
     * reflects the transformation operations performed
     * on the Processiong.js transformation functions.
     * Sane people will most likely invoke this after a
     * resetMatrix().
     * The Interpose parameter is the Interpose library that
     * knows how to capture calls to Processing.js.  Get that
     * library at 
     * https://www.khanacademy.org/cs/i/5992045511049216
     */
    RMatrix2D.mirrorPJS = function(Interpose) {
        var mirror = RMatrix2D.new();
        
        /*
         * This interpose-function invokes the equivalent
         * method in the RMatrix2D object "mirror" that 
         * reflects the Processing.js rendering matrix.
         */
        var mirrorMatrixOp = function() {
            mirror[this.name].apply(mirror, arguments);
        };
        
        /* Interpose plumbing... */
        for (var f = RMatrix2D.transforms, i = f.length; --i >= 0; ) {
            Interpose.push(f[i], mirrorMatrixOp);
        }
        
        /* Undo method.  The method itself disappears... */
        mirror.breakMirror = function() {
            /* unplumb Interpose... */
            for (var f = RMatrix2D.transforms, i = f.length; --i >= 0; ) {
                Interpose.pop(f[i]);
            }
            delete mirror.breakMirror;
        };
        
        return mirror;
    };
    
    return RMatrix2D; // last edited 2019.02.23
})();  /* RMatrix2D is a constructor */
var mops, rm = new RMatrix2D();
var tuX = new PVector();  /* transformed unit X vector */
var tuY = new PVector();  /* transformed unit Y vector */
var tX = width/2;
var tY = height/2;
var sideLen = min(0.55*width, 0.55*height);
sideLen = 2 * round(sideLen / 2);
var halfLen = sideLen / 2;
var deepestZ = -halfLen;
var nodes=[]; //Blank Node Set for working with
var faces=[];
var moving=false;
var sideColor,tails,heads; //Coin Dependencies.
var coinRadius; //Please Change for Development
var face_round = function(img,w,h) {
    ellipse(0, 0, w, h);
    if (img) {
        image(img, 0, 0);
    }
}; //Face for Round Circle
for (var i = 0; i < faces.length; i++) {
    var osb = faces[i].screen = createGraphics(sideLen, sideLen, JAVA2D);
    osb.background(faces[i].color);
    osb.textSize(0.7 * sideLen);
    osb.fill(0, 0, 0);
    osb.textAlign(CENTER, CENTER);
    osb.text(faces[i].label, (sideLen-1)/2, (sideLen-1)/2);
} //Face Code (Cube)
var drawFace = function(f, drawEdges) {
    var n = f.nodes, ox = n[0].x, oy = n[0].y;
    rm.set( (n[1].x - ox) / sideLen, (n[3].x - ox) / sideLen, ox,
            (n[1].y - oy) / sideLen, (n[3].y - oy) / sideLen, oy);
        
    pushMatrix();
    RMatrix2D.applyMatrix(rm);
    image(f.screen, 0, 0);
    popMatrix();
    
    for (var j = n.length-1, i = 0; drawEdges && i < n.length; j = i, i++) {
        line(n[j].x, n[j].y, n[i].x, n[i].y);
    }
}; //Cubes!
var rotate3D = function(theta, xProp, yProp) {
    var sinTheta = sin(theta);
    var cosTheta = cos(theta);

    for (var minZ = Infinity, n = 0; n < nodes.length; n++) {
        var node = nodes[n];
        var x = node[xProp];
        var y = node[yProp];
        node[xProp] = x * cosTheta - y * sinTheta;
        node[yProp] = y * cosTheta + x * sinTheta;
        minZ = min(minZ, node.z);
    }
    
    return minZ;
}; //Rotation Class

// Rotate shape around the z-axis.
var rotateZ3D = function(theta) {
    return rotate3D(theta, "x", "y");
};

// Rotate shape around the y-axis & find minimum z value.
var rotateY3D = function(theta) {
    return rotate3D(theta, "z", "x");
};

// Rotate shape around the x-axis & find minimum z value.
var rotateX3D = function(theta) {
    return rotate3D(theta, "y", "z");
};

var faceIsVisible = function(face) {
    for (var nodes = face.nodes, n = 0; n < nodes.length; n++) {
        if (nodes[n].z === deepestZ) {
            return false;
        }
    }
    return true;
};

var drawCylinderSide = function(tail) {
    var origin = tail.nodes[0], ox = origin.x, oy = origin.y;
    var heads = tail.other.nodes[0], hx = heads.x, hy = heads.y;
    /* 
     * Set vector (dx, dy) to the perpendicular of the vector 
     * between heads and origin, with magnitude coinRadius.
     */
    var dy = hx - ox;
    var dx = oy - hy;
    var mag = sqrt(dx*dx + dy*dy);
    dx *= coinRadius/mag;
    dy *= coinRadius/mag;
    /* The view of the side is a quad with the following coordinates: */
    var ax = ox + dx, ay = oy + dy;
    var bx = hx + dx, by = hy + dy;
    var cx = hx - dx, cy = hy - dy;
    var dx = ox - dx, dy = oy - dy;
    
    /* Fill(?) with noStroke... */
    pushStyle();
    noStroke();
    quad(ax, ay, bx, by, cx, cy, dx, dy);
    popStyle();
    
    /* Stroke (maybe)... */
    line(ax, ay, bx, by);
    line(cx, cy, dx, dy);
};

var drawCylinder = function() {
    var bottom = (heads.nodes[0].z > tails.nodes[0].z) ? tails : heads;
    var top  = (bottom === tails) ? heads : tails;
    
    var n = top.nodes, origin = n[0];
    /* tuX = normailze(n[1] - n[0]) */
    tuX.set(n[1]);
    tuX.sub(origin);
    tuX.div(coinRadius);
    /* tuY = normalize(n[2] - n[0]) */
    tuY.set(n[2]);
    tuY.sub(origin);
    tuY.div(coinRadius);
    rm.set( tuX.x, tuY.x, 0,
            tuX.y, tuY.y, 0 );
    mops = rm.reconstitute(mops);
    
    fill(sideColor);
    drawFace(bottom);
    drawCylinderSide(bottom);
    fill(top.color);
    drawFace(top, top.img);
};
