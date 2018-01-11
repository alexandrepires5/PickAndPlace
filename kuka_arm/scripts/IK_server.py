#!/usr/bin/env python

# Copyright (C) 2017 Electric Movement Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import radians
from sympy import *
import numpy as np


def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
		print ("Starting program")
    	q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8') #this is theta
    	d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
    	a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
    	alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')
		#DH parameters definition
    	dh = {alpha0: 0, a0: 0, d1: 0.75, 
    	alpha1: -pi/2, a1: 0.35, d2: 0, q2: q2 - pi/2,
    	alpha2: 0, a2: 1.25, d3: 0,
    	alpha3: -pi/2, a3: -0.054, d4: 1.5,
    	alpha4: pi/2, a4: 0, d5: 0,
    	alpha5: -pi/2, a5: 0, d6: 0,
    	alpha6: 0, a6: 0, d7: 0.303, q7: 0 } 
		xc, yc, zc = symbols('xc yc zc')
		A, B, C = symbols('A B C')
		yaw, pitch, roll = symbols('yaw pitch roll')
		t1, t2, t3, t4, t5, t6 = symbols('t1 t2 t3 t4 t5 t6')
		r11, r12, r13, r21, r22, r23, r31, r32, r33 = symbols('r11 r12 r13 r21 r22 r23 r31 r32 r33')
		#thetas calculation
		t1 = atan2(yc, xc)
		t2 = pi/2 - atan2(zc, yc) - acos((C*C + B*B - A*A) / (2*B*C))
		t3 = pi/2 - (acos((A*A + C*C - B*B) / (2*A*C)) + 0.036)
		t4 = atan2(r33, -r13)
		t5 = atan2(sqrt(r13 * r13 + r33 * r33), r23)
		t6 = atan2(-r22, r21)

		print ("Finish variables init")
    	R_z = Matrix([[cos(np.pi), -sin(np.pi), 0, 0],
                  [sin(np.pi), cos(np.pi), 0, 0],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1]])
    
    
    	R_y = Matrix([[cos(-np.pi / 2), 0, sin(-np.pi / 2), 0],
                  [0, 1, 0, 0],
                  [-sin(-np.pi / 2), 0, cos(-np.pi / 2), 0],
                  [0, 0, 0, 1]])
		R_corr = R_z * R_y
        print ("Correction matrix calculated")
    	#------------------------------------------------
    
    	#Calculate nx, ny and nz parameters to calculate WC
    	R_z_yaw = Matrix([[cos(yaw), -sin(yaw), 0, 0],
                  [sin(yaw), cos(yaw), 0, 0],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1]])
    
    	R_y_pitch = Matrix([[cos(pitch), 0, sin(pitch), 0],
                  [0, 1, 0, 0],
                  [-sin(pitch), 0, cos(pitch), 0],
                  [0, 0, 0, 1]])
    
    	R_x_roll = Matrix([[1, 0, 0, 0],
                  [0, cos(roll), -sin(roll), 0],
                  [0, sin(roll), cos(roll), 0],
                  [0, 0, 0, 1]])

		T0_1 = Matrix([[cos(q1), -sin(q1), 0, a0],
		      [sin(q1)*cos(alpha0), cos(q1)*cos(alpha0), -sin(alpha0), -sin(alpha0)*d1],
		      [sin(q1)*sin(alpha0), cos(q1)*sin(alpha0), cos(alpha0), cos(alpha0)*d1],
		      [0, 0, 0, 1]])

    	T0_1 = T0_1.subs(dh)


    	T1_2 = Matrix([[cos(q2), -sin(q2), 0, a1],
		      [sin(q2)*cos(alpha1), cos(q2)*cos(alpha1), -sin(alpha1), -sin(alpha1)*d2],
		      [sin(q2)*sin(alpha1), cos(q2)*sin(alpha1), cos(alpha1), cos(alpha1)*d2],
		      [0, 0, 0, 1]])
    	T1_2 = T1_2.subs(dh)

    	T2_3 = Matrix([[cos(q3), -sin(q3), 0, a2],
		      [sin(q3)*cos(alpha2), cos(q3)*cos(alpha2), -sin(alpha2), -sin(alpha2)*d3],
		      [sin(q3)*sin(alpha2), cos(q3)*sin(alpha2), cos(alpha2), cos(alpha2)*d3],
		      [0, 0, 0, 1]])
    	T2_3 = T2_3.subs(dh)

    	T3_4 = Matrix([[cos(q4), -sin(q4), 0, a3],
		      [sin(q4)*cos(alpha3), cos(q4)*cos(alpha3), -sin(alpha3), -sin(alpha3)*d4],
		      [sin(q4)*sin(alpha3), cos(q4)*sin(alpha3), cos(alpha3), cos(alpha3)*d4],
		      [0, 0, 0, 1]])

    	T3_4 = T3_4.subs(dh)


    	T4_5 = Matrix([[cos(q5), -sin(q5), 0, a4],
		      [sin(q5)*cos(alpha4), cos(q5)*cos(alpha4), -sin(alpha4), -sin(alpha4)*d5],
		      [sin(q5)*sin(alpha4), cos(q5)*sin(alpha4), cos(alpha4), cos(alpha4)*d5],
		      [0, 0, 0, 1]])
    	T4_5 = T4_5.subs(dh)

    	T5_6 = Matrix([[cos(q6), -sin(q6), 0, a5],
		      [sin(q6)*cos(alpha5), cos(q6)*cos(alpha5), -sin(alpha5), -sin(alpha5)*d6],
		      [sin(q6)*sin(alpha5), cos(q6)*sin(alpha5), cos(alpha5), cos(alpha5)*d6],
		      [0, 0, 0, 1]])
    	T5_6 = T5_6.subs(dh)

    	T6_G = Matrix([[cos(q7), -sin(q7), 0, a6],
		      [sin(q7)*cos(alpha6), cos(q7)*cos(alpha6), -sin(alpha6), -sin(alpha6)*d7],
		      [sin(q7)*sin(alpha6), cos(q7)*sin(alpha6), cos(alpha6), cos(alpha6)*d7],
		      [0, 0, 0, 1]])

    	T6_G = T6_G.subs(dh)	
		print ("Transformation matrices defined")
		R0_3 = T0_1[0:3, 0:3] * T1_2[0:3, 0:3] * T2_3[0:3, 0:3]
        Rrpy = R_z_yaw * R_y_pitch * R_x_roll * R_corr

		T0_2 = T0_1 * T1_2
        T0_3 = T0_2 * T2_3
        T0_4 = T0_3 * T3_4
    	T0_5 = T0_4 * T4_5
    	T0_6 = T0_5 * T5_6
    	T0_G = T0_6 * T6_G

    	T_total = T0_G * R_corr
        # Initialize service response
        joint_trajectory_list = []
        end_effector_error = []
        for x in xrange(0, len(req.poses)):
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()
	    	print ("IK starting")
	    	# Extract end-effector position and orientation from request
	    	# px,py,pz = end-effector position
	    	# roll, pitch, yaw = end-effector orientation
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z
            #Orientations for end effector obtained with euler_from_quaternion transformation
            (r, p, y) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, req.poses[x].orientation.y,
                    req.poses[x].orientation.z, req.poses[x].orientation.w])
            #evaluate roll pitch yaw (r p y) angles to obtain rotation matrix from end effector with the correction aswell
            Rrpy = R_z_yaw.evalf(subs={yaw: y}) * R_y_pitch.evalf(subs={pitch: p}) * R_x_roll.evalf(subs={roll: r}) * R_corr
            print ("Orientation obtained")
	    	[nx, ny, nz] = [Rrpy[0,2], Rrpy[1,2], Rrpy[2,2]]
            
            l = 0.303
    	    d6 = 0
    
    	    WC_x = (px - (d6 + l)*nx)
    	    WC_y = (py - (d6 + l)*ny)
    	    WC_z = (pz - (d6 + l)*nz)	
            print ("Wrist position calculated")
            #sides a b c obtained from DH parameters + extra calculation from b side, since we need to consider it is a 3-D work environment.
	    	a = 1.501
	    	b = sqrt(((sqrt(WC_x**2 + WC_y**2) - 0.35))**2 + pow((WC_z - 0.75),2))
    	    c = 1.25

	    	#theta 1 to theta 3 calculation
	    	theta1 = t1.evalf(subs={xc: WC_x, yc: WC_y})
    	    theta2 = t2.evalf(subs={zc: WC_z - 0.75, yc: sqrt(WC_x**2 + WC_y**2)-0.35, A: a, B: b, C: c})
    	    theta3 = t3.evalf(subs={A: a, B: b, C: c})
	    
            #theta 4 to 6 calculation from rotation matrix
	    	R_with_thetas = R0_3.evalf(subs={q1: theta1, q2:theta2, q3:theta3})
            R3_6 = R_with_thetas.inv("LU")* Rrpy[0:3, 0:3]

            theta4 = t4.evalf(subs={r13: R3_6[0,2], r33: R3_6[2,2]})
    	    theta5 = t5.evalf(subs={r13: R3_6[0,2], r33: R3_6[2,2], r23: R3_6[1,2]})
            theta6 = t6.evalf(subs={r22: R3_6[1,1], r21: R3_6[1,0]})
            #calculation of end effector position from input angles
            T_final = T_total.evalf(subs={q1: theta1, q2:theta2, q3:theta3, q4: theta4, q5: theta5, q6: theta6})
            [ee_x, ee_y, ee_z] = [T_final[0,3], T_final[1,3], T_final[2,3]]
            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
	    	joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
	    	joint_trajectory_list.append(joint_trajectory_point)
            #calculate end effector error
            ee_x_e = abs(ee_x-req.poses[x].position.x)
            ee_y_e = abs(ee_y-req.poses[x].position.y)
            ee_z_e = abs(ee_z-req.poses[x].position.z)
            ee_offset = sqrt(ee_x_e**2 + ee_y_e**2 + ee_z_e**2)
            end_effector_error.append(ee_offset)

        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        print ("Return IK response")
        print ("End effector error for all points of trajectory is:", end_effector_error)
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()