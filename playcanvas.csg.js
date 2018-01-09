var  ii = 0.5;
pc.Vec3.TransformCoordinates = function(e,i){
    var r = pc.Vec3.ZERO.clone();
    return pc.Vec3.TransformCoordinatesToRef(e,i,r),r;
};

pc.Vec3.TransformCoordinatesToRef = function(e, t, i) {
    var r = e.x * t.data[0] + e.y * t.data[4] + e.z * t.data[8] + t.data[12]
      , n = e.x * t.data[1] + e.y * t.data[5] + e.z * t.data[9] + t.data[13]
      , o = e.x * t.data[2] + e.y * t.data[6] + e.z * t.data[10] + t.data[14]
      , s = e.x * t.data[3] + e.y * t.data[7] + e.z * t.data[11] + t.data[15];

    i.x = r / s,
    i.y = n / s,
    i.z = o / s
};


 pc.Vec3.TransformNormal = function(e, i) {
     var r = pc.Vec3.ZERO.clone();
     
     return pc.Vec3.TransformNormalToRef(e, i, r),
            r;
 };
        
pc.Vec3.TransformNormalToRef = function(e, t, i) {
    var r = e.x * t.data[0] + e.y * t.data[4] + e.z * t.data[8]
      , n = e.x * t.data[1] + e.y * t.data[5] + e.z * t.data[9]
      , o = e.x * t.data[2] + e.y * t.data[6] + e.z * t.data[10];
    i.x = r,
    i.y = n,
    i.z = o
};

var currentCSGMeshId = 0;
var attributeMap = { position: pc.SEMANTIC_POSITION, normal: pc.SEMANTIC_NORMAL, tangent: pc.SEMANTIC_TANGENT, blendWeight: pc.SEMANTIC_BLENDWEIGHT, blendIndices: pc.SEMANTIC_BLENDINDICES, color: pc.SEMANTIC_COLOR, texCoord0: pc.SEMANTIC_TEXCOORD0, texCoord1: pc.SEMANTIC_TEXCOORD1, texCoord2: pc.SEMANTIC_TEXCOORD2, texCoord3: pc.SEMANTIC_TEXCOORD3, texCoord4: pc.SEMANTIC_TEXCOORD4, texCoord5: pc.SEMANTIC_TEXCOORD5, texCoord6: pc.SEMANTIC_TEXCOORD6, texCoord7: pc.SEMANTIC_TEXCOORD7 };
var Vertex = (function() {

    function Vertex(pos, normal, uv) {
        this.pos = pos;
        this.normal = normal;
        this.uv = uv;
    };

    Vertex.prototype.clone = function() {
        return new Vertex(this.pos.clone(), this.normal.clone(), this.uv.clone());
    };

    Vertex.prototype.flip = function() {
        this.normal = this.normal.scale(-1); ///////// not found;
    };

    Vertex.prototype.interpolate = function(other, t) {

        var lerpPos = this.pos.clone();
        lerpPos.lerp(this.pos, other.pos, t);

        var lerpNor = this.normal.clone();
        lerpNor.lerp(this.normal, other.normal, t);

        var lerpUv = this.uv.clone();
        lerpUv.lerp(this.uv, other.uv, t);

        var lerpb;
        return new Vertex(lerpPos, lerpNor, lerpUv);
    };

    return Vertex;
}());


var Plane = (function() {
    function Plane(normal, w) {
        this.normal = normal;
        this.w = w;
    };

    Plane.FromPoints = function(a, b, c) {
        var cc = c.clone(),
            bb = b.clone(),
            aa = a.clone();
        
        var v0 = cc.sub(aa);
        var v1 = bb.sub(aa);

        if (v0.lengthSq() === 0 || v1.lengthSq() === 0)
            return null;
        var n = new pc.Vec3();
        n.cross(v0, v1);
        n.normalize();

        return new Plane(n, n.dot(a));
    };

    Plane.prototype.clone = function() {
        return new Plane(this.normal.clone(), this.w);
    };

    Plane.prototype.flip = function() {
        this.normal.scale(-1);
        this.w = -this.w;
    };
    Plane.prototype.splitPolygon = function(polygon, coplanarFront, coplanarBack, front, back) {
        var COPLANAR = 0;
        var FRONT = 1;
        var BACK = 2;
        var SPANNING = 3;

        var polygonType = 0;
        var types = new Array();
        var t,dot,w = this.w;
        
        for(var i = 0,ii=polygon.vertices.length; i<ii ;i++){
            dot = this.normal.dot(polygon.vertices[i].pos);
            t = dot -w;
            var type = (t < -Plane.EPSILON) ? BACK : (t > Plane.EPSILON) ? FRONT : COPLANAR;
            polygonType |= type;
            types.push(type);
        }

        switch(polygonType) {
            case COPLANAR:
                var nor = this.normal.clone();
                (nor.dot(polygon.plane.normal) > 0 ? coplanarFront : coplanarBack).push(polygon);
                break;
            case FRONT:
                front.push(polygon);
                break;
            case BACK:
                back.push(polygon);
                break;
            case SPANNING:
                var f = [],
                    b = [];
                
                for(var i = 0 , ii=polygon.vertices.length; i< ii;i++){
                    
                    var j = (i + 1) % polygon.vertices.length;

                    var ti = types[i],
                        tj = types[j];

                    var vi = polygon.vertices[i],
                        vj = polygon.vertices[j];

                    if (ti !== BACK)
                        f.push(vi);

                    if (ti !== FRONT)
                        b.push(ti !== BACK ? vi.clone() : vi);

                    if ((ti | tj) === SPANNING) {
                        
                        var vjPos = new pc.Vec3();
                        t = (this.w - this.normal.dot(vi.pos)) / this.normal.dot(vjPos.sub2(vj.pos,vi.pos));
                        var v = vi.interpolate(vj, t);
                        f.push(v);
                        b.push(v.clone());
                    }
                }

                var poly;
                if (f.length >= 3) {
                    poly = new Polygon(f, polygon.shared);
                    if (poly.plane)
                        front.push(poly);
                }
                if (b.length >= 3) {
                    poly = new Polygon(b, polygon.shared);
                    if (poly.plane)
                        back.push(poly);
                }
                break;
        }
    };
    Plane.EPSILON = 1e-5 ;
    
    return Plane;
}());


var Polygon = (function() {
    function Polygon(vertices, shared) {
        this.vertices = vertices;
        this.shared = shared;
        this.plane = Plane.FromPoints(vertices[0].pos, vertices[1].pos, vertices[2].pos);
    }

    Polygon.prototype.clone = function() {
        var vertices = this.vertices.map(function(v) { return v.clone(); });
        return new Polygon(vertices, this.shared);
    };

    Polygon.prototype.flip = function() {
        this.vertices.reverse().map(function(v) { v.flip(); });
        this.plane.flip();
    };

    return Polygon;
}());



var Node = (function() {
    function Node(polygons) {
        this.plane = null;
        this.front = null;
        this.back = null;
        this.polygons = [];

        if (polygons) {
            this.build(polygons);
        }
    };

    Node.prototype.clone = function() {
        var node = new Node();
        node.plane = this.plane && this.plane.clone();
        node.front = this.front && this.front.clone();
        node.back = this.back && this.back.clone();
        node.polygons = this.polygons.map(function(p) { return p.clone(); });
        return node;
    };

    Node.prototype.invert = function() {
        for (var i = 0; i < this.polygons.length; i++) {
            this.polygons[i].flip();
        }
        if (this.plane) {
            this.plane.flip();
        }
        if (this.front) {
            this.front.invert();
        }
        if (this.back) {
            this.back.invert();
        }
        var temp = this.front;
        this.front = this.back;
        this.back = temp;
    };

    // Recursively remove all polygons in `polygons` that are inside this BSP
    // tree.
    Node.prototype.clipPolygons = function(polygons) {
        if (!this.plane)
            return polygons.slice();
        var front = [],
            back = [];
        for (var i = 0,ii= polygons.length; i <ii; i++) {
            this.plane.splitPolygon(polygons[i], front, back, front, back);
        }
        if (this.front) {
            front = this.front.clipPolygons(front);
        }
        if (this.back) {
            back = this.back.clipPolygons(back);
        } else {
            back = [];
        }
        return front.concat(back);
    };

    Node.prototype.clipTo = function(bsp) {
        this.polygons = bsp.clipPolygons(this.polygons);
        if (this.front) {
            this.front.clipTo(bsp);
        }
        if (this.back)
            this.back.clipTo(bsp);
    };
    Node.prototype.allPolygons = function() {
        var polygons = this.polygons.slice();

        if (this.front)
            polygons = polygons.concat(this.front.allPolygons());

        if (this.back)
            polygons = polygons.concat(this.back.allPolygons());

        return polygons;
    };
    // Build a BSP tree out of `polygons`. When called on an existing tree, the
    // new polygons are filtered down to the bottom of the tree and become new
    // nodes there. Each set of polygons is partitioned using the first polygon
    // (no heuristic is used to pick a good split).
    Node.prototype.build = function(polygons) {
        if (!polygons.length)
            return;
        if (!this.plane)
            this.plane = polygons[0].plane.clone();
        var front = [],
            back = [];
        for (var i = 0,ii=polygons.length; i < ii; i++) {
            this.plane.splitPolygon(polygons[i], this.polygons, this.polygons, front, back);
        }
        if (front.length) {
            if (!this.front)
                this.front = new Node();
            this.front.build(front);
        }
        if (back.length) {
            if (!this.back)
                this.back = new Node();
            this.back.build(back);
        }
    };

    return Node;
}());


var CSG = (function() {
    function CSG() {
        this.polygons = [];
    };


    CSG.FromMesh = function(entity) {
        var  normal, uv, position, polygon, polygons = new Array(),
            vertices=[];
        var matrix, meshPosition, meshRotation, meshRotationQuaternion = null,
            meshScaling,scalePercent;
        
        var meshInstances = entity.model.model.meshInstances;
        var meshInstance= meshInstances[0],
            mesh = meshInstance.mesh,
            graphNode = entity.model.model.graph,
            indexBuffer = mesh.indexBuffer[0],
            buffer = mesh.vertexBuffer,
            iterator = new pc.VertexIterator(buffer),
            format = buffer.format;
        
                
        var positions = [],
            normals = [],
            uvs = [];
        
        var attribute, attributeName;

        if (graphNode instanceof pc.GraphNode) {
            
            matrix = graphNode.getWorldTransform();
            meshPosition = graphNode.position.clone();
            meshRotation = graphNode.localEulerAngles.clone();
            if (graphNode.rotation) {
                meshRotationQuaternion = graphNode.localRotation.clone();
            }
            meshScaling = graphNode.localScale.clone();
            
        } else
            throw 'Playcanvas.CSG: Wrong Mesh type, must be Playcanvas.Mesh';

        var numVertices = iterator.element[attributeMap.position].array.length /  format.size*4;
        
        for (var j = 0,jj= numVertices; j < jj; j++) {

            for (var k = 0; k < format.elements.length; k++) {
                attributeName = format.elements[k].name;
                var index = iterator.element[attributeName].index;
                var array = iterator.element[attributeName].array;

                switch (attributeName) {
                    case attributeMap.position:
                        positions.push(array[index]);
                        positions.push(array[index + 1]);
                        positions.push(array[index + 2]);
                        break;

                    case attributeMap.normal:
                        normals.push(array[index]);
                        normals.push(array[index + 1]);
                        normals.push(array[index + 2]);
                        break;

                    case attributeMap.texCoord0:
                        uvs.push(array[index]);
                        uvs.push(array[index + 1]);
                        break;
                }

            }
            iterator.next();
        }
        iterator.end();


        var indices= new Uint16Array(indexBuffer.storage, 0);

        for (var j = 0,jj = indices.length; j < jj; j +=3) {
            vertices = [];

            for(var n = 0,nn=3;n<nn;n++){
                var sourcePos = new pc.Vec3(positions[indices[j+n] * 3], positions[indices[j+n] * 3 + 1], positions[indices[j+n] * 3 + 2]);
                var sourceNor = new pc.Vec3(normals[indices[j+n]*3], normals[indices[j+n] * 3 + 1], normals[indices[j+n] * 3 + 2]);
                var sourceUv = new pc.Vec2(uvs[indices[j+n] * 2], uvs[indices[j+n] * 2 + 1]);
                var pos =pc.Vec3.TransformCoordinates(sourcePos,matrix);
                var nor =pc.Vec3.TransformNormal(sourceNor,matrix);

                vertices.push(new Vertex(pos, nor, sourceUv));
            }

            polygon = new Polygon(vertices);

            if (polygon.plane)
                polygons.push(polygon);
        }
        
        var csg = CSG.FromPolygons(polygons);
        csg.matrix = matrix;
        csg.position = meshPosition;
        csg.rotation = meshRotation;
        csg.scaling = meshScaling;
        csg.rotationQuaternion = meshRotationQuaternion;
        csg.entity = entity;
        currentCSGMeshId++;
        return csg;
    };

    CSG.FromPolygons = function(polygons) {
        var csg = new CSG();
        csg.polygons = polygons;
        return csg;
    };

    CSG.prototype = {
        clone: function() {
            var csg = new CSG();
            csg.polygons = this.polygons.map(function(p) { return p.clone(); });
            csg.copyTransformAttributes(this);
            return csg;
        },
        intersect: function(csg) {
            var a = new Node(this.clone().polygons);
            var b = new Node(csg.clone().polygons);
            a.clipTo(b);
            b.clipTo(a);
            b.invert();
            b.clipTo(a);
            b.invert();
            a.build(b.allPolygons());
            return CSG.FromPolygons(a.allPolygons()).copyTransformAttributes(this);
        },
        intersectInPlace: function(csg) {
            var a = new Node(this.polygons);
            var b = new Node(csg.polygons);
            a.clipTo(b);
            b.clipTo(a);
            b.invert();
            b.clipTo(a);
            b.invert();
            a.build(b.allPolygons());
            this.polygons = a.allPolygons();
        },
        subtract: function(csg) {
            var b = new Node(this.clone().polygons);
            var a = new Node(csg.clone().polygons);
            a.invert();
            a.clipTo(b);
            b.clipTo(a);
            b.invert();
            b.clipTo(a);
            b.invert();
            a.build(b.allPolygons());
            // a.invert(); 
            return CSG.FromPolygons(a.allPolygons()).copyTransformAttributes(this);
        },
        subtractInPlace: function(csg) {
            var a = new Node(this.polygons);
            var b = new Node(csg.polygons);
            a.invert();
            a.clipTo(b);
            b.clipTo(a);
            b.invert();
            b.clipTo(a);
            b.invert();
            a.build(b.allPolygons());
            a.invert();
            this.polygons = a.allPolygons();
        },

        union: function(csg) {
            var a = new Node(this.clone().polygons);
            var b = new Node(csg.clone().polygons);
            a.invert();
            b.clipTo(a);
            b.invert();
            a.clipTo(b);
            b.clipTo(a);
            a.build(b.allPolygons());
            a.invert();
            return CSG.FromPolygons(a.allPolygons()).copyTransformAttributes(this);
        },
        unionInPlace: function(csg) {
            var a = new Node(this.polygons);
            var b = new Node(csg.polygons);
            a.invert();
            b.clipTo(a);
            b.invert();
            a.clipTo(b);
            b.clipTo(a);
            a.build(b.allPolygons());
            a.invert();
            this.polygons = a.allPolygons();
        },
        // Return a new BABYLON.CSG solid with solid and empty space switched. This solid is
        // not modified.
        inverse: function() {
            var csg = this.clone();
            csg.inverseInPlace();
            return csg;
        },
        inverseInPlace: function() {
            this.polygons.map(function(p) { p.flip(); });
        },
        copyTransformAttributes: function(csg) {
            this.matrix = csg.matrix;
            this.position = csg.position;
            this.rotation = csg.rotation;
            this.scaling = csg.scaling;
            this.rotationQuaternion = csg.rotationQuaternion;
            this.entity = csg.entity;
            return this;
        },

        buildMeshGeometry: function() {
            var matrix = this.matrix.clone();
            matrix.invert();
            var mesh,
                vertices = [],
                indices = [],
                normals = [],
                uvs = [],
                vertex = pc.Vec3.ZERO,
                normal = pc.Vec3.ZERO,
                uv = pc.Vec2.ZERO,
                polygons = this.polygons,
                polygonIndices = [0, 0, 0],
                polygon, vertice_dict = {},
                vertex_idx, currentIndex = 0,
                subMesh_dict = {},
                subMesh_obj;
            
            
            for (var i = 0, il = polygons.length; i < il; i++) {
                polygon = polygons[i];
                // Building SubMeshes

                for (var j = 2, jl = polygon.vertices.length; j < jl; j++) {
                    polygonIndices[0] = 0;
                    polygonIndices[1] = j - 1;
                    polygonIndices[2] = j;
                    for (var k = 0; k < 3; k++) {
                        vertex = polygon.vertices[polygonIndices[k]].pos.clone();
                        normal = polygon.vertices[polygonIndices[k]].normal.clone();
                        
                        uv = polygon.vertices[polygonIndices[k]].uv.clone();

                        var localVertex = pc.Vec3.TransformCoordinates(vertex,matrix);
                        var localNormal = pc.Vec3.TransformNormal(normal,matrix);

                        vertex_idx = vertice_dict[localVertex.x + ',' + localVertex.y + ',' + localVertex.z];
                        // Check if 2 points can be merged
                        if (!(typeof vertex_idx !== 'undefined' &&
                                normals[vertex_idx * 3] === localNormal.x &&
                                normals[vertex_idx * 3 + 1] === localNormal.y &&
                                normals[vertex_idx * 3 + 2] === localNormal.z &&
                                uvs[vertex_idx * 2] === uv.x &&
                                uvs[vertex_idx * 2 + 1] === uv.y)) {
                            
                            // localVertex.scale(100);
                            // normal.scale(100);
                            vertices.push(localVertex.x, localVertex.y, localVertex.z);
                            normals.push(normal.x, normal.y, normal.z);
                            uvs.push(uv.x, uv.y);
                            
                            vertex_idx = vertice_dict[localVertex.x + ',' + localVertex.y + ',' + localVertex.z] = (vertices.length / 3) - 1;
                        }
                        indices.push(vertex_idx);
                        currentIndex++;
                    }
                }
            }
            
            
            var options = {
                normals : normals,
                uvs     : uvs,
                indices : indices
            };
            
            mesh = pc.createMesh(pc.app.graphicsDevice, vertices, options);
          
            return mesh;
        },

        toMesh: function() {
            var mesh = this.buildMeshGeometry();

            var node = new pc.GraphNode();
            node.setLocalEulerAngles(this.rotation);
            node.setLocalRotation(this.rotationQuaternion);
            node.setLocalScale(this.scaling);
            node.setLocalPosition(this.position);
            
            var material = new pc.StandardMaterial();
            var meshInstance = new pc.MeshInstance(node,mesh,material);
            
            
            var model = new pc.Model();
            model.graph = node; 
            model.meshInstances.push(meshInstance);
            
            
            
            model.generateWireframe();
            
            meshInstance.renderStyle = pc.RENDERSTYLE_WIREFRAME;
            pc.app.scene.removeModel(this.entity.model.model);
            this.entity.model.model.destroy();
            
            this.entity.removeComponent('model');

            var e = this.entity;
            e.addComponent('model');
            e.model.model = model;
            
        }
    };
    return CSG;
}());