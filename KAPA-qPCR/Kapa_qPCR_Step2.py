#!/usr/bin/env python
"""
Step 2 of KAPA qPCR

Original Metadata:
__author__ = "Joe Brown, Greg Sjogren, Diane Luo"
__credits__ = ["Joe Brown", "Greg Sjogren", "Diane Luo"]
__version__ = "1.0"
__date__ = "2021-09-02"
"""

import math
import json
from opentrons import protocol_api, types #needed for trying to set specific pipette movements before pipetting
from opentrons import simulate


def get_values(*names):
    _all_values = json.loads("""{"temp_deck":"temperature module gen2", 
    "pipette_type":"p20_multi_gen2",
    "pipette_mount":"right",
    "pipette_type_2":"p300_multi_gen2",
    "pipette_mount_2":"left",
    "sample_number":9,
    "sample_volume":4,
    "master_mix_volume":6.2,
    "set_temperature":4}""")
    return [_all_values[n] for n in names]

metadata = {
    'protocolName': 'Kapa_Illumina Library qPCR',
    'author': 'Greg Sjogren, Joe Brown, Diane Luo',
    'source': 'Kapa KR0405 v9.17 Protocol',
    'apiLevel': '2.10', 'softwareLevel': '4.50'
    }


def run(protocol_context):

    [temp_deck, pipette_type, pipette_mount, pipette_type_2, pipette_mount_2, sample_number, sample_volume,
     master_mix_volume, set_temperature] = get_values(  # noqa: F821
        "temp_deck", "pipette_type", "pipette_mount", "pipette_type_2", "pipette_mount_2", "sample_number",
        "sample_volume", "master_mix_volume", "set_temperature"
    )
    # Set Speed of Z Axis
    protocol_context.max_speeds['Z'] = 100
    
    # Set tip touch-off to true. Needs to be turned on at start of any script where tiptouching is used
    protocol_context.touch_tip = True
     
    # Load Temperature deck in Slot 10
    temp_deck = protocol_context.load_module(temp_deck, '10')
   
   # Populate temp_deck with 96 well PCR Reagent Plate
    temp_plate = temp_deck.load_labware(
        'biorad_96_wellplate_200ul_pcr')
        
    # Load 384 well plate in Slot 1  
    qPCR_plate = protocol_context.load_labware(
        'biorad_384_wellplate_50ul', '1', 'qPCR plate')
        
    # Load Custom Labware BioRad 384 well plate (Jupyter Simulations Only). Note that the .json file needs to br uploaded into the Jupyter Notebook main directory. There is conflciting info that the .json file needs to be in a "labware" folder on there; this is tested and untrue. 
    #with open('biorad_384_wellplate_50ul.json') as labware_file:
        #labware_def = json.load(labware_file)
    #reagent_container = protocol_context.load_labware_from_definition(labware_def, '1')
        
    # Load 1-10k 96 well dilution plate in Slot 2 
    dilution_10k_plate = protocol_context.load_labware(
        'biorad_96_wellplate_200ul_pcr', '2', 'dilution 10k plate')
        
    # Load 1-20k 96 well dilution plate in Slot 3 
    dilution_20k_plate = protocol_context.load_labware(
        'biorad_96_wellplate_200ul_pcr', '3', 'dilution 20k plate')
        
    #Math to make loops work for samples variables
    col_num = math.ceil(sample_number/8)# IE the total # columns you will be processing. 
    col_num_384_10k = math.ceil(sample_number/8*2)#Setting the endpoint for the slices below, based on the # of sample columns you need to process for the 10K Dilution Plate
    col_num_384_20k = math.ceil(sample_number/8*2)+12#Setting the endpoint for the slices below, based on the # of sample columns you need to process for the 20K Dilution Plate
    #col_offset = (168 - sample_number) / 8# Leave in if we change to dipsnese standards dynamically. 
    #standard_col_num = math.ceil(sample_number/8+col_offset)#Leave in if we change to dipsnese standards dynamically.
    output_qPCR_quad_1_10k = [col for col in qPCR_plate.rows()[0][:col_num_384_10k:2]]#Start at Row A, 1st column (leaving no text before the first: means it will start in the first column; this could start with 0. Pipette SampleNumber/8 (rounded up) # of Columns, skipping a column in between each Sample Column that is pipetted.  
    output_qPCR_quad_2_10k = [col for col in qPCR_plate.rows()[0][1:col_num_384_10k:2]]#Start at Row A, 2nd column (1 is used since Python is 0 based). Pipette SampleNumber/8 (rounded up) # of Columns, skipping a column in between each Sample Column that is pipetted.
    output_qPCR_quad_3_10k = [col for col in qPCR_plate.rows()[1][:col_num_384_10k:2]]#Start at Row B, 1st column (leaving no text before the first: means it will start in the first column; this could start with 0. Pipette SampleNumber/8 (rounded up) # of Columns, skipping a column in between each Sample Column that is pipetted.
    output_qPCR_quad_1_20k = [col for col in qPCR_plate.rows()[0][12:col_num_384_20k:2]]#Start at Row A, 13th column (12 is used since Python is 0 based). Pipette SampleNumber/8 (rounded up) # of Columns, skipping a column in between each Sample Column that is pipetted.
    output_qPCR_quad_2_20k = [col for col in qPCR_plate.rows()[0][13:col_num_384_20k:2]]##Start at Row A, 14th column (13 is used since Python is 0 based). Pipette SampleNumber/8 (rounded up) # of Columns, skipping a column in between each Sample Column that is pipetted.
    output_qPCR_quad_3_20k = [col for col in qPCR_plate.rows()[1][12:col_num_384_20k:2]]#Start at Row B, 13th column (12 is used since Python is 0 based). Pipette SampleNumber/8 (rounded up) # of Columns, skipping a column in between each Sample Column that is pipetted.
    output_standards_quad_4 = [col for col in qPCR_plate.rows()[1][1:6:2]]#Start at Row B, 2nd column (1 is used since Python is 0 based), and proceed untilyou reach column 6. In this case, MM and standards go into B2, B4, and B6. 
    
    # Use only 20 uL tips per sample in this protocol, tip rack goes in slots 4-9 and 11
    total_tips = (((col_num*8)*12)+48)# Max Number of Tips for 48 sample run. 
    tiprack_num = math.ceil(total_tips/96)
    slots = ['4', '5', '6', '7', '8', '9', '11'][:tiprack_num]

    tip_name = 'opentrons_96_tiprack_20ul'
    tipracks = [protocol_context.load_labware(tip_name, slot) for slot in slots]

    # Telling Pippete Mount (right_pipette, in this case 20 ul multichannel) to use 20 uL tips
    pipette = protocol_context.load_instrument(
        pipette_type, pipette_mount, tip_racks=tipracks)
     
    

    # Define Reagent Source Columns
    master_mix_col_1 = temp_plate.wells()[0]#IE in Column 1
    master_mix_col_2 = temp_plate.wells()[8] # IE in Column 2
    standards_col_4 = temp_plate.wells()[32]# IE in **Column 5**
    
    #Math for User Deck Preparation
    load_tips = ((col_num*8)*12)+48
    load_tip_boxes = math.ceil(load_tips/96)
    col_1_MM = math.ceil((((col_num*8)*3)*master_mix_volume)/8+20)
    col_2_MM = math.ceil(((((col_num*8)*3)+24)*master_mix_volume)/8+20)
    col_4_STDs = (sample_volume*3)+20
    
    protocol_context.set_rail_lights(True)
    #############################################################################################################################################################################
    ##########################################################User Deck Preparation Prompts######################################################################################
    #############################################################################################################################################################################
    protocol_context.pause("If not already present, load Temperature Module (Gen2) onto deck grid 10, plug in power & USB and press ON button. Wipe Down Deck with 70% Ethanol")
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("""Load {} 20 uL Tip Boxes onto deck positions in the following order: 4,5,6,7,8,9,11 """.format(str(load_tip_boxes)))
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("Please load BioRad Hardshell 384-well qPCR Plate onto deck grid 1.Tape down with lab tape so side touches do not lift plate off of th deck.")
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("Load BioRad Hardshell 96-well 10k Dilution Plate onto deck grid 2.")
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("Load BioRad Hardshell 96-well 20k Dilution Plate onto deck grid 3.")
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("""Set up BioRad Hardshell 96-well Reagent Plate as follows: Column 1- pipette {} uL Master Mix into all column wells.""".format(str(col_1_MM)))
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("""Set up BioRad Hardshell 96-well Reagent Plate as follows: Column 2- pipette {} uL Master Mix into all column wells.""".format(str(col_2_MM)))
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("""Set up BioRad Hardshell 96-well Reagent Plate as follows: Column 5- pipette {} uL Standards & NTCs into all column wells.""".format(str(col_4_STDs)))
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("Load BioRad Hardshell 96-well Reagent Plate Plate onto the Temperature Module (Gen2) on deck grid 10.")
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("Review Deck Layout Photo. Make Sure All Plates are unsealed and tip rack overs are removed. Once you click resume, pipetting will begin!")
    ##############################################################################################################################################################################
    ##############################################################################################################################################################################
    
    #Turn off Deck Lights as MasterMix is Light Sensitive
    protocol_context.set_rail_lights(False)
    
    #Set Temperature of Temperature Module to User Defined Variable
    temp_deck.set_temperature(set_temperature)
    
   #Define Master Mix Final Aspirate Volume
    master_mix_volume = (master_mix_volume)-0.3# subtraction of 0.8 is to account for pipette overdelivering an average of 0.35 uL per trasnsfer
    
    # Define Master Mix Flow Rates
    pipette.flow_rate.aspirate = 6.2
    pipette.flow_rate.dispense = 6.2
    pipette.flow_rate.blow_out = 20
    
    # Define Aspiration Position
    if col_num == 6:
        master_mix_Well_Bottom=7
    if col_num == 5:
        master_mix_Well_Bottom=6.4
    if col_num == 4:
        master_mix_Well_Bottom=5.8
    if col_num == 3:
        master_mix_Well_Bottom=5.2
    if col_num == 2:
        master_mix_Well_Bottom=4.6
    if col_num == 1:
        master_mix_Well_Bottom=4

    # Dispense qPCR MM for Dilution Plate 1 Rep 1
    for target in output_qPCR_quad_1_10k:
        pipette.pick_up_tip()
        pipette.well_bottom_clearance.aspirate = master_mix_Well_Bottom
        pipette.well_bottom_clearance.dispense = 2
        pipette.aspirate(master_mix_volume, master_mix_col_1)
        pipette.dispense(master_mix_volume, target)
        protocol_context.delay(seconds=3)
        pipette.blow_out(target.bottom(2.5))
        pipette.touch_tip(v_offset=-5)
        pipette.drop_tip()
        master_mix_Well_Bottom=master_mix_Well_Bottom-0.4
        
     # Define Aspiration Position
    if col_num == 6:
        master_mix_Well_Bottom=4.6
    if col_num == 5:
        master_mix_Well_Bottom=4.4
    if col_num == 4:
        master_mix_Well_Bottom=4.2
    if col_num == 3:
        master_mix_Well_Bottom=4
    if col_num == 2:
        master_mix_Well_Bottom=3.8
    if col_num == 1:
        master_mix_Well_Bottom=3.6
    
    # Dispense qPCR MM for Dilution Plate 1 Rep 2
    for target in output_qPCR_quad_2_10k:
        pipette.pick_up_tip()
        pipette.well_bottom_clearance.aspirate = master_mix_Well_Bottom
        pipette.well_bottom_clearance.dispense = 2
        pipette.aspirate(master_mix_volume, master_mix_col_1)
        pipette.dispense(master_mix_volume, target)
        protocol_context.delay(seconds=3)
        pipette.blow_out(target.bottom(2.5))
        pipette.touch_tip(v_offset=-5)
        pipette.drop_tip()
        master_mix_Well_Bottom=master_mix_Well_Bottom-0.4
        
    # Define Aspiration Position
    if col_num == 6:
        master_mix_Well_Bottom=2.2
    if col_num == 5:
        master_mix_Well_Bottom=2.4
    if col_num == 4:
        master_mix_Well_Bottom=2.6
    if col_num == 3:
        master_mix_Well_Bottom=2.8
    if col_num == 2:
        master_mix_Well_Bottom=3.0
    if col_num == 1:
        master_mix_Well_Bottom=3.2
        
    # Dispense qPCR MM for Dilution Plate 1 Rep 3
    for target in output_qPCR_quad_3_10k:
        pipette.pick_up_tip()
        pipette.well_bottom_clearance.aspirate = master_mix_Well_Bottom
        pipette.well_bottom_clearance.dispense = 2
        pipette.aspirate(master_mix_volume, master_mix_col_1)
        pipette.dispense(master_mix_volume, target)
        protocol_context.delay(seconds=3)
        pipette.blow_out(target.bottom(2.5))
        pipette.touch_tip(v_offset=-5)
        pipette.drop_tip()
        master_mix_Well_Bottom=master_mix_Well_Bottom-0.4  
    
    # Define Aspiration Position for MM Column 2
    if col_num == 6:
        master_mix_Well_Bottom=7.4
    if col_num == 5:
        master_mix_Well_Bottom=6.8
    if col_num == 4:
        master_mix_Well_Bottom=6.2
    if col_num == 3:
        master_mix_Well_Bottom=5.6
    if col_num == 2:
        master_mix_Well_Bottom=5
    if col_num == 1:
        master_mix_Well_Bottom=4.4
        
    # Dispense qPCR MM for Dilution Plate 2 Rep 1
    for target in output_qPCR_quad_1_20k:
        pipette.pick_up_tip()
        pipette.well_bottom_clearance.aspirate = master_mix_Well_Bottom
        pipette.well_bottom_clearance.dispense = 2
        pipette.aspirate(master_mix_volume, master_mix_col_2)
        pipette.dispense(master_mix_volume, target)
        protocol_context.delay(seconds=3)
        pipette.blow_out(target.bottom(2.5))
        pipette.touch_tip(v_offset=-5)
        pipette.drop_tip()
        master_mix_Well_Bottom=master_mix_Well_Bottom-0.4
        
     # Define Aspiration Position
    if col_num == 6:
        master_mix_Well_Bottom=5
    if col_num == 5:
        master_mix_Well_Bottom=4.8
    if col_num == 4:
        master_mix_Well_Bottom=4.6
    if col_num == 3:
        master_mix_Well_Bottom=4.4
    if col_num == 2:
        master_mix_Well_Bottom=4.2
    if col_num == 1:
        master_mix_Well_Bottom=4
        
    # Dispense qPCR MM for Dilution Plate 2 Rep 2
    for target in output_qPCR_quad_2_20k:
        pipette.pick_up_tip()
        pipette.well_bottom_clearance.aspirate = master_mix_Well_Bottom
        pipette.well_bottom_clearance.dispense = 2
        pipette.aspirate(master_mix_volume, master_mix_col_2)
        pipette.dispense(master_mix_volume, target)
        protocol_context.delay(seconds=3)
        pipette.blow_out(target.bottom(2.5))
        pipette.touch_tip(v_offset=-5)
        pipette.drop_tip()
        master_mix_Well_Bottom=master_mix_Well_Bottom-0.4
        
    # Define Aspiration Position
    if col_num == 6:
        master_mix_Well_Bottom=2.6
    if col_num == 5:
        master_mix_Well_Bottom=2.8
    if col_num == 4:
        master_mix_Well_Bottom=3
    if col_num == 3:
        master_mix_Well_Bottom=3.2
    if col_num == 2:
        master_mix_Well_Bottom=3.4
    if col_num == 1:
        master_mix_Well_Bottom=3.6
        
    # Dispense qPCR MM for Dilution Plate 2 Rep 3
    for target in output_qPCR_quad_3_20k:
        pipette.pick_up_tip()
        pipette.well_bottom_clearance.aspirate = master_mix_Well_Bottom
        pipette.well_bottom_clearance.dispense = 2
        pipette.aspirate(master_mix_volume, master_mix_col_2)
        pipette.dispense(master_mix_volume, target)
        protocol_context.delay(seconds=3)
        pipette.blow_out(target.bottom(2.5))
        pipette.touch_tip(v_offset=-5)
        pipette.drop_tip()
        master_mix_Well_Bottom=master_mix_Well_Bottom-0.4
        
     # Define Aspiration Position
    master_mix_Well_Bottom=1
    
    # Dispense qPCR MM for Standards & NTCs
    for target in output_standards_quad_4:
        pipette.pick_up_tip()
        pipette.well_bottom_clearance.aspirate = master_mix_Well_Bottom
        pipette.well_bottom_clearance.dispense = 2
        pipette.aspirate(master_mix_volume, master_mix_col_2)
        pipette.dispense(master_mix_volume, target)
        protocol_context.delay(seconds=3)
        pipette.blow_out(target.bottom(2.5))
        pipette.touch_tip(v_offset=-5)
        pipette.drop_tip()
    
    protocol_context.set_rail_lights(True)
    #############################################################################################################################################################################
    ##########################################################Prompt to Centrifuge qPCR BioRad Hard Shell 384-well plate#########################################################
    #############################################################################################################################################################################  
    protocol_context.pause("Remove qPCR BioRad Hard Shell 384-well plate from deck grid 1. Seal Plate.")
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("Centrifuge plate briefly to remove all bubbles. ReLoad qPCR BioRad Hard Shell 384-well plate onto deck grid 1. Remove seal. Re-Tape Plate to Deck. Empty Trash! Once you click resume, pipetting will begin!")
    protocol_context.set_rail_lights(False)
    ##############################################################################################################################################################################
    ##############################################################################################################################################################################
    
    protocol_context.set_rail_lights(False)
    
    # Define Sample Transfer Flow Rates
    pipette.flow_rate.aspirate = 4
    pipette.flow_rate.dispense = 4
    pipette.flow_rate.blow_out = 20
    
    #Define Sample Final Aspirate Volume
    sample_volume = (sample_volume)+0.1 # additon of 0.1 is to account for pipette underdelivering an average of 0.125 uL per transfer
    
    # Define Dilution Plate Aspiration Position Replicate 1
    dilution_plate_aspirate_Height=5
    
    # Transfer Dilution Plate 10K to Quadrant 1 of qPCR Plate. This Code structure will allow you to aspirate from 1 variable source column into 1 variable destination column. This is controlled by the col_num counter. 
    for e in range(col_num):
        if e == 0:
            dilution_10K_source = dilution_10k_plate['A1']
            qPCR_Dest_10K_Quad_1 = qPCR_plate['A1']
        if e == 1:
            dilution_10K_source = dilution_10k_plate['A2']
            qPCR_Dest_10K_Quad_1 = qPCR_plate['A3']
        if e == 2:
            dilution_10K_source = dilution_10k_plate['A3']
            qPCR_Dest_10K_Quad_1 = qPCR_plate['A5']
        if e == 3:
            dilution_10K_source = dilution_10k_plate['A4']
            qPCR_Dest_10K_Quad_1 = qPCR_plate['A7']
        if e == 4:
            dilution_10K_source = dilution_10k_plate['A5']
            qPCR_Dest_10K_Quad_1 = qPCR_plate['A9']
        if e == 5:
            dilution_10K_source = dilution_10k_plate['A6']
            qPCR_Dest_10K_Quad_1 = qPCR_plate['A11']
        
        pipette.pick_up_tip()
        pipette.well_bottom_clearance.aspirate = dilution_plate_aspirate_Height
        pipette.well_bottom_clearance.dispense = 2
        pipette.aspirate(sample_volume, dilution_10K_source)
        pipette.dispense(sample_volume, qPCR_Dest_10K_Quad_1)
        pipette.well_bottom_clearance.aspirate = 2
        pipette.mix(3, sample_volume, qPCR_Dest_10K_Quad_1)
        protocol_context.delay(seconds=3)
        pipette.blow_out(qPCR_Dest_10K_Quad_1.bottom(3.5))
        pipette.touch_tip(v_offset=-5)
        pipette.drop_tip()
        
    # Define Dilution Plate Aspiration Position Replicate 2
    dilution_plate_aspirate_Height=4.6
    
      # Transfer Dilution Plate 10K to Quadrant 2 of qPCR Plate
    for e in range(col_num):
        if e == 0:
            dilution_10K_source = dilution_10k_plate['A1']
            qPCR_Dest_10K_Quad_2 = qPCR_plate['A2']
        if e == 1:
            dilution_10K_source = dilution_10k_plate['A2']
            qPCR_Dest_10K_Quad_2 = qPCR_plate['A4']
        if e == 2:
            dilution_10K_source = dilution_10k_plate['A3']
            qPCR_Dest_10K_Quad_2 = qPCR_plate['A6']
        if e == 3:
            dilution_10K_source = dilution_10k_plate['A4']
            qPCR_Dest_10K_Quad_2 = qPCR_plate['A8']
        if e == 4:
            dilution_10K_source = dilution_10k_plate['A5']
            qPCR_Dest_10K_Quad_2 = qPCR_plate['A10']
        if e == 5:
            dilution_10K_source = dilution_10k_plate['A6']
            qPCR_Dest_10K_Quad_2 = qPCR_plate['A12']
        
        pipette.pick_up_tip()
        pipette.well_bottom_clearance.aspirate = dilution_plate_aspirate_Height
        pipette.well_bottom_clearance.dispense = 2
        pipette.aspirate(sample_volume, dilution_10K_source)
        pipette.dispense(sample_volume, qPCR_Dest_10K_Quad_2)
        pipette.well_bottom_clearance.aspirate = 2
        pipette.mix(3, sample_volume, qPCR_Dest_10K_Quad_2)
        protocol_context.delay(seconds=3)
        pipette.blow_out(qPCR_Dest_10K_Quad_2.bottom(3.5))
        pipette.touch_tip(v_offset=-5)
        pipette.drop_tip()
        
        
   # Define Dilution Plate Aspiration Position Replicate 3  
    dilution_plate_aspirate_Height=4.2
    
      # Transfer Dilution Plate 10K to Quadrant 3 of qPCR Plate
    for e in range(col_num):
        if e == 0:
            dilution_10K_source = dilution_10k_plate['A1']
            qPCR_Dest_10K_Quad_3 = qPCR_plate['B1']
        if e == 1:
            dilution_10K_source = dilution_10k_plate['A2']
            qPCR_Dest_10K_Quad_3 = qPCR_plate['B3']
        if e == 2:
            dilution_10K_source = dilution_10k_plate['A3']
            qPCR_Dest_10K_Quad_3 = qPCR_plate['B5']
        if e == 3:
            dilution_10K_source = dilution_10k_plate['A4']
            qPCR_Dest_10K_Quad_3 = qPCR_plate['B7']
        if e == 4:
            dilution_10K_source = dilution_10k_plate['A5']
            qPCR_Dest_10K_Quad_3 = qPCR_plate['B9']
        if e == 5:
            dilution_10K_source = dilution_10k_plate['A6']
            qPCR_Dest_10K_Quad_3 = qPCR_plate['B11']
        
        pipette.pick_up_tip()
        pipette.well_bottom_clearance.aspirate = dilution_plate_aspirate_Height
        pipette.well_bottom_clearance.dispense = 2
        pipette.aspirate(sample_volume, dilution_10K_source)
        pipette.dispense(sample_volume, qPCR_Dest_10K_Quad_3)
        pipette.well_bottom_clearance.aspirate = 2
        pipette.mix(3, sample_volume, qPCR_Dest_10K_Quad_3)
        protocol_context.delay(seconds=3)
        pipette.blow_out(qPCR_Dest_10K_Quad_3.bottom(3.5))
        pipette.touch_tip(v_offset=-5)
        pipette.drop_tip()
             
        # Define Dilution Plate Aspiration Position Replicate 1
    dilution_plate_aspirate_Height=7
        
    # Transfer Dilution Plate 20K to Quadrant 1 of qPCR Plate
    for e in range(col_num):
        if e == 0:
            dilution_20K_source = dilution_20k_plate['A1']
            qPCR_Dest_20K_Quad_1 = qPCR_plate['A13']
        if e == 1:
            dilution_20K_source = dilution_20k_plate['A2']
            qPCR_Dest_20K_Quad_1 = qPCR_plate['A15']
        if e == 2:
            dilution_20K_source = dilution_20k_plate['A3']
            qPCR_Dest_20K_Quad_1 = qPCR_plate['A17']
        if e == 3:
            dilution_20K_source = dilution_20k_plate['A4']
            qPCR_Dest_20K_Quad_1 = qPCR_plate['A19']
        if e == 4:
            dilution_20K_source = dilution_20k_plate['A5']
            qPCR_Dest_20K_Quad_1 = qPCR_plate['A21']
        if e == 5:
            dilution_20K_source = dilution_20k_plate['A6']
            qPCR_Dest_20K_Quad_1 = qPCR_plate['A23']
        
        pipette.pick_up_tip()
        pipette.well_bottom_clearance.aspirate = dilution_plate_aspirate_Height
        pipette.well_bottom_clearance.dispense = 2
        pipette.aspirate(sample_volume, dilution_20K_source)
        pipette.well_bottom_clearance.aspirate = 2
        pipette.dispense(sample_volume, qPCR_Dest_20K_Quad_1)
        pipette.mix(3, sample_volume, qPCR_Dest_20K_Quad_1)
        protocol_context.delay(seconds=3)
        pipette.blow_out(qPCR_Dest_20K_Quad_1.bottom(3.5))
        pipette.touch_tip(v_offset=-5)
        pipette.drop_tip()
        
    # Define Dilution Plate Aspiration Position Replicate 2
    dilution_plate_aspirate_Height=6.6
        
      # Transfer Dilution Plate 20K to Quadrant 2 of qPCR Plate
    for e in range(col_num):
        if e == 0:
            dilution_20K_source = dilution_20k_plate['A1']
            qPCR_Dest_20K_Quad_2 = qPCR_plate['A14']
        if e == 1:
            dilution_20K_source = dilution_20k_plate['A2']
            qPCR_Dest_20K_Quad_2 = qPCR_plate['A16']
        if e == 2:
            dilution_20K_source = dilution_20k_plate['A3']
            qPCR_Dest_20K_Quad_2 = qPCR_plate['A18']
        if e == 3:
            dilution_20K_source = dilution_20k_plate['A4']
            qPCR_Dest_20K_Quad_2 = qPCR_plate['A20']
        if e == 4:
            dilution_20K_source = dilution_20k_plate['A5']
            qPCR_Dest_20K_Quad_2 = qPCR_plate['A22']
        if e == 5:
            dilution_20K_source = dilution_20k_plate['A6']
            qPCR_Dest_20K_Quad_2 = qPCR_plate['A24']
        
        pipette.pick_up_tip()
        pipette.well_bottom_clearance.aspirate = dilution_plate_aspirate_Height
        pipette.well_bottom_clearance.dispense = 2
        pipette.aspirate(sample_volume, dilution_20K_source)
        pipette.dispense(sample_volume, qPCR_Dest_20K_Quad_2)
        pipette.well_bottom_clearance.aspirate = 2
        pipette.mix(3, sample_volume, qPCR_Dest_20K_Quad_2)
        protocol_context.delay(seconds=3)
        pipette.blow_out(qPCR_Dest_20K_Quad_2.bottom(3.5))
        pipette.touch_tip(v_offset=-5)
        pipette.drop_tip()
        
        
    # Define Dilution Plate Aspiration Position Replicate 3
    dilution_plate_aspirate_Height=6.2
    
      # Transfer Dilution Plate 10K to Quadrant 3 of qPCR Plate
    for e in range(col_num):
        if e == 0:
            dilution_20K_source = dilution_20k_plate['A1']
            qPCR_Dest_20K_Quad_3 = qPCR_plate['B13']
        if e == 1:
            dilution_20K_source = dilution_20k_plate['A2']
            qPCR_Dest_20K_Quad_3 = qPCR_plate['B15']
        if e == 2:
            dilution_20K_source = dilution_20k_plate['A3']
            qPCR_Dest_20K_Quad_3 = qPCR_plate['B17']
        if e == 3:
            dilution_20K_source = dilution_20k_plate['A4']
            qPCR_Dest_20K_Quad_3 = qPCR_plate['B19']
        if e == 4:
            dilution_20K_source = dilution_20k_plate['A5']
            qPCR_Dest_20K_Quad_3 = qPCR_plate['B21']
        if e == 5:
            dilution_20K_source = dilution_20k_plate['A6']
            qPCR_Dest_20K_Quad_3 = qPCR_plate['B23']
        
        pipette.pick_up_tip()
        pipette.well_bottom_clearance.aspirate = dilution_plate_aspirate_Height
        pipette.well_bottom_clearance.dispense = 2
        pipette.aspirate(sample_volume, dilution_20K_source)
        pipette.dispense(sample_volume, qPCR_Dest_20K_Quad_3)
        pipette.well_bottom_clearance.aspirate = 2
        pipette.mix(3, sample_volume, qPCR_Dest_20K_Quad_3)
        protocol_context.delay(seconds=3)
        pipette.blow_out(qPCR_Dest_20K_Quad_3.bottom(3.5))
        pipette.touch_tip(v_offset=-5)
        pipette.drop_tip()
        
    # Define Dilution Plate Aspiration Position for Standards & NTCs
    dilution_plate_aspirate_Height=2.1
    sample_volume = (sample_volume-0.2) # to account for 0.2 uL overdispense based on standards being kept @4C due to low concetration.
        
     # Dispense Standards & NTCs into qPCR Plate
    for target in output_standards_quad_4:
        pipette.pick_up_tip()
        pipette.well_bottom_clearance.aspirate = dilution_plate_aspirate_Height
        pipette.well_bottom_clearance.dispense = 2
        pipette.aspirate(sample_volume, standards_col_4)
        pipette.dispense(sample_volume, target)
        pipette.well_bottom_clearance.aspirate = 2
        pipette.mix(3, sample_volume, target)
        protocol_context.delay(seconds=3)
        pipette.blow_out(target.bottom(3.5))
        pipette.touch_tip(v_offset=-5)
        pipette.drop_tip()

    protocol_context.set_rail_lights(True)
        
    # Turn off Temperature Deck
    temp_deck.deactivate()
    
    ################################################################################################################################################################################
    ######################################################################Deck Unloading Instructions###############################################################################
    ################################################################################################################################################################################
    protocol_context.pause("Please remove Tip Waste from deck grid 12 to biohazard bin.")
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("Remove Dilution Plates from deck grid 2 & 3 to biohazard bin.")
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("Remove Reagent Plate from deck grid 10 Temperature Module to biohazard bin.")
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("Seal qPCR plate with MicroAmp Optical Adhesive Cover and remove from deck.")
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("Remove Tip Boxes from deck. Wipe Down Deck with 70% Ethanol. Once you click 'Resume', pipette head will raise and script will end!")
    ################################################################################################################################################################################
    ################################################################################################################################################################################
    
    #protocol_context.set_rail_lights(False) #if this is left on, then on your next run, you cannot keep light on even if you turn it on via code. The light will flash but not stay on.
    protocol_context.set_rail_lights(False) # is above note true?? Test with real samples. Test with Closing opentrons App, then turning back on. If that doesnt work, try hard booting the Robot, and see if it works.    
    
if __name__ == '__main__':
    protocol = simulate.get_protocol_api('2.10')
    run(protocol)
    for line in protocol.commands():
        print(line)
