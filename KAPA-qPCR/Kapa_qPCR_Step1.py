#!/usr/bin/env python
"""
Step 1 of KAPA qPCR

Original Metadata:
__author__ = "Joe Brown, Greg Sjogren, Diane Luo"
__credits__ = ["Joe Brown", "Greg Sjogren", "Diane Luo"]
__version__ = "1.0"
__date__ = "2021-08-23"
"""


import math
import json
from opentrons import protocol_api, types #needed for trying to set specific pipette movements before pipetting
from opentrons import simulate


def get_values(*names):
    _all_values = json.loads("""{"temp_deck":"temperature module gen2", 
    "pipette_type":"p300_multi_gen2",
    "pipette_mount":"left",
    "pipette_type_2":"p20_multi_gen2",
    "pipette_mount_2":"right",
    "sample_number":9,
    "starting_sample_volume":12,
    "dilution_volume_1":998,
    "sample_volume_1":2,
    "dilution_volume_2":95,
    "sample_volume_2":5,
    "dilution_volume_3":40,
    "sample_volume_3":40,
    "set_temperature":4}""")
    return [_all_values[n] for n in names]

metadata = {
    'protocolName': 'Kapa_Illumina Library qPCR Step 1',
    'author': 'Greg Sjogren, Joe Brown, Diane Luo',
    'source': 'Kapa KR0405 v9.17 Protocol',
    'apiLevel': '2.10', 'softwareLevel': '4.50'
    }

def run(protocol_context):

    [temp_deck, pipette_type, pipette_mount, pipette_type_2, pipette_mount_2, sample_number, starting_sample_volume, sample_volume_1, dilution_volume_1, sample_volume_2, dilution_volume_2, sample_volume_3, dilution_volume_3, set_temperature] = get_values(  # noqa: F821
        "temp_deck", "pipette_type", "pipette_mount", "pipette_type_2", "pipette_mount_2", "sample_number", "starting_sample_volume",
        "sample_volume_1", "dilution_volume_1", "sample_volume_2", "dilution_volume_2", "sample_volume_3", "dilution_volume_3", "set_temperature"
    )
    # Set Speed of Z Axis
    protocol_context.max_speeds['Z'] = 100
    
    # Check and see If Lights are On; Turn Lights On if Currently Off; Need to Troubleshoot this
    #protocol_context.rail_lights_on
    #original_light_status = protocol_context.rail_lights_on
    #if original_light_status == True:
        #protocol_context.comment("Lights are already on!")
    #if original_light_status == False:
        #protocol_context.comment("Turning Lights on!")
        #protocol_context.set_rail_lights(True)
    #else:
        #protocol_context.comment("You have defied boolean logic!")
        
    # Turn Lights On
    protocol_context.set_rail_lights(True)
        
    # Load Agilent 4 well 73 mL Reagent Reservoir
    reagent_container = protocol_context.load_labware(
        'agilent_4_well_73_ml_reagent_reservoir', '4')
        
    # Load Agilent 4 well 73 mL Reagent Reservoir (Jupyter Simulations ONly). Note that the .json file needs to br uploaded into the Jupyter Notebook main directory. There is conflciting info that the .json file needs to be in a "labware" folder on there; this is tested and untrue. 
    #with open('agilent_4_well_73_ml_reagent_reservoir.json') as labware_file:
        #labware_def = json.load(labware_file)
    #reagent_container = protocol_context.load_labware_from_definition(labware_def, '4')
     
    # Load Temperature deck in Slot 10
    temp_deck = protocol_context.load_module(temp_deck, '10')
   
   # Populate temp_deck with final 1-20K 96 well PCR Reagent Plate
    dilution_20k_plate = temp_deck.load_labware(
        'biorad_96_wellplate_200ul_pcr')
        
    # Load Inital Sample Plate 96 well plate in Slot 1  
    sample_plate = protocol_context.load_labware(
        'biorad_96_wellplate_200ul_pcr', '1', 'Sample plate')
        
    # Load 1-500 2 mL 96 well dilution plate in Slot 2 
    dilution_500_plate = protocol_context.load_labware(
        'perkinelmer_96_wellplate_2000ul', '2', 'dilution plate 1')
        
    # Load 1-10k 96 well dilution plate in Slot 3 
    dilution_10k_plate = protocol_context.load_labware(
        'biorad_96_wellplate_200ul_pcr', '3', 'dilution plate 2 10k plate')
        
     # Define samples variables
    col_num = math.ceil(sample_number/8)# IE the total # columns you will be processing. 
    samples = [col for col in sample_plate.rows()[0][:col_num]]
    samples_dilution_1 = [col for col in dilution_500_plate.rows()[0][:col_num]]
    samples_dilution_2 = [col for col in dilution_10k_plate.rows()[0][:col_num]]
    samples_dilution_3 = [col for col in dilution_20k_plate.rows()[0][:col_num]]
    samples_10k = [col for col in dilution_10k_plate.rows()[0][:col_num]]
    samples_20k = [col for col in dilution_20k_plate.rows()[0][:col_num]]
        
    # Use four 300 uL tips per sample, tip racks go in slots 5-8 300uL tips
    total_tips = (col_num*8)*4
    tiprack_num = math.ceil(total_tips/96)
    slots = ['5', '6', '7', '8'][:tiprack_num]
    
    tip_name = 'opentrons_96_tiprack_300ul'
    tipracks = [protocol_context.load_labware(tip_name, slot) for slot in slots]

    # Use one 20 uL tips per sample, tip rack goes in slot 10 20 uL tips
    total_tips_2 = (col_num*8)*2
    tiprack_num_2 = math.ceil(total_tips_2/96)
    slots_2 = ['9'][:tiprack_num_2]

    tip_name_2 = 'opentrons_96_tiprack_20ul'
    tipracks_2 = [protocol_context.load_labware(tip_name_2, slot_2) for slot_2 in slots_2]
    
    # Telling Pippete Mount (left_pipette, in this case 300 ul multichannel) to use 300 uL tips
    pipette = protocol_context.load_instrument(
        pipette_type, pipette_mount, tip_racks=tipracks)
        
    # Telling Pippete Mount 2(right pipette, in this case 20 ul multichannel) to use 20 uL tips 
    pipette_2 = protocol_context.load_instrument(
        pipette_type_2, pipette_mount_2, tip_racks=tipracks_2)
     
    # Define reagents and liquid waste
    dilution_buffer = reagent_container.wells()[0]
    
    # Starting Sample Volume Error Handling
    if starting_sample_volume == 12:
        starting_aspirate_height = 0.2
        pipette_2.well_bottom_clearance.aspirate = starting_aspirate_height
    else: 
        protocol_context.pause("Starting sample volume is out of range, use 12uL only. Please re-run protocol script and correct volume.")
    
    #Define User Deck Preparation
    load_tips = ((col_num*8)*4)
    load_300_tip_boxes = math.ceil(load_tips/96)
    load_tips2 = ((col_num*8)*2)
    load_20_tip_boxes = math.ceil(load_tips2/96)
    col_1_Dilution_Buffer = math.ceil((((col_num*8)*(dilution_volume_1+dilution_volume_2+dilution_volume_3)+5000)/1000))
    
    #############################################################################################################################################################################
    ##########################################################User Deck Preparation Prompts######################################################################################
    #############################################################################################################################################################################
    protocol_context.pause("If not already present, load Temperature Module (Gen2) onto deck grid 10, plug in power & USB and press ON button. Wipe Down Deck with 70% Ethanol.")
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("Label and Load BioRad 96-well Hard Shell 20k Dilution Plate onto Temperature Module on deck grid 10")
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("Load  1 20 uL Tip Box onto deck grid 9")
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("""Load {} 300 uL Tip Boxes onto deck positions in the following order: 5,6,7,8 """.format(str(load_300_tip_boxes)))
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("""1. Load Agilent 73 mL Reagent Reservoir onto deck grid 4.    2. Pipette {} mL Dilution Buffer (10mM TrisHCL, 0.5% Tween20) into Well A1.""".format(str(col_1_Dilution_Buffer)))
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("Label BioRad 96-well Hard Shell Initial Sample Plate (Initial Dilution).Vortex Plate for 1 minute at Speed 10. Cnetirufe for 500 x g for 2 minutes.Load BioRad 96-well Hard Shell Initial Sample Plate onto deck grid 1.")
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("Label and Load PE Pipetting Microplate, 2mL DW SQ 96-well plate onto deck grid 2")
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("Label and Load BioRad 96-well Hard Shell 10k Dilution Plate onto deck grid 3")
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("Review Deck Layout Photo. Make Sure All Plates are unsealed and tip rack overs are removed. Once you click resume, pipetting will begin!")
  
    ##############################################################################################################################################################################
    ##############################################################################################################################################################################
    
    #Turn off Deck Lights 
    protocol_context.set_rail_lights(False)
    
    # Define Well Bottom for Reagent Reservoir (agilent_4_well_73_ml_reagent_reservoir) to dilution_500_plate (perkinelmer_96_wellplate_2000ul) 
    if col_num == 6:
        Dilution_Buffer_Well_Bottom=25
    if col_num == 5:
        Dilution_Buffer_Well_Bottom=21
    if col_num == 4:
        Dilution_Buffer_Well_Bottom=17
    if col_num == 3:
        Dilution_Buffer_Well_Bottom=13
    if col_num == 2:
        Dilution_Buffer_Well_Bottom=9
    if col_num == 1:
        Dilution_Buffer_Well_Bottom=5
    
   #Define Dilution Final Aspirate Volume beyond tip capacity
    aspirate_volume = (dilution_volume_1/4)-0.8# subtraction of 0.8 is to account for pipette overdelivering an average of 0.8 uL per trasnsfer
    
    # Change Flow Rates
    pipette.flow_rate.aspirate = 94 #speeding up so dispense can be very low for accuracy. 94 choosen as its pipette default speed.  
    pipette.flow_rate.dispense = 22.5 #250 and 95 leave a decent amount of volume in tips when blowout. at 250 the volume is moving too fast for the buffer to escape surface tension. Still happens at 30, but less so. 15 looks great , but is SLOWWWWW. take the middle between the two, which is 22.5. doesnt look as good, but takes 4 minutes per column to fill up. Fairly repoducbile too, +/- 2 uL, and at best Opentrons promises +/-1.5 uL 
    pipette.flow_rate.blow_out = 299

    # Dispense Dilution Buffer to PE Pipetting Microplate 2mL DW SQ 96-well plate
    for target in samples_dilution_1:
        protocol_context.max_speeds['Z'] = 200 # Set Speed of Z Axis
        pipette.pick_up_tip()
        pipette.well_bottom_clearance.aspirate = Dilution_Buffer_Well_Bottom
        pipette.well_bottom_clearance.dispense = Dilution_Buffer_Well_Bottom
        pipette.mix(1, aspirate_volume, dilution_buffer)
        pipette.well_bottom_clearance.aspirate = Dilution_Buffer_Well_Bottom
        pipette.well_bottom_clearance.dispense = 4
        protocol_context.max_speeds['Z'] = 20 # Set Speed of Z Axis
        pipette.aspirate(aspirate_volume, dilution_buffer)
        protocol_context.max_speeds['Z'] = 200 # Set Speed of Z Axis
        pipette.dispense(aspirate_volume, target)
        protocol_context.delay(seconds=3)
        pipette.blow_out(target.bottom(4.5))
        pipette.touch_tip(v_offset=-5)
        pipette.well_bottom_clearance.aspirate = Dilution_Buffer_Well_Bottom
        pipette.well_bottom_clearance.dispense = Dilution_Buffer_Well_Bottom
        pipette.mix(1, aspirate_volume, dilution_buffer)
        pipette.well_bottom_clearance.aspirate = Dilution_Buffer_Well_Bottom
        pipette.well_bottom_clearance.dispense = 8
        protocol_context.max_speeds['Z'] = 20 # Set Speed of Z Axis
        pipette.aspirate(aspirate_volume, dilution_buffer)
        protocol_context.max_speeds['Z'] = 200 # Set Speed of Z Axis
        pipette.dispense(aspirate_volume, target)
        protocol_context.delay(seconds=3)
        pipette.blow_out(target.bottom(8.5))
        pipette.touch_tip(v_offset=-5)
        pipette.well_bottom_clearance.aspirate = Dilution_Buffer_Well_Bottom
        pipette.well_bottom_clearance.dispense = Dilution_Buffer_Well_Bottom
        pipette.mix(1, aspirate_volume, dilution_buffer)
        pipette.well_bottom_clearance.aspirate = Dilution_Buffer_Well_Bottom
        pipette.well_bottom_clearance.dispense = 12
        protocol_context.max_speeds['Z'] = 20 # Set Speed of Z Axis
        pipette.aspirate(aspirate_volume, dilution_buffer)
        protocol_context.max_speeds['Z'] = 200 # Set Speed of Z Axis
        pipette.dispense(aspirate_volume, target)
        protocol_context.delay(seconds=3)
        pipette.blow_out(target.bottom(12.5))
        pipette.touch_tip(v_offset=-5)
        pipette.well_bottom_clearance.aspirate = Dilution_Buffer_Well_Bottom
        pipette.well_bottom_clearance.dispense = Dilution_Buffer_Well_Bottom
        pipette.mix(1, aspirate_volume, dilution_buffer)
        pipette.well_bottom_clearance.aspirate = Dilution_Buffer_Well_Bottom
        pipette.well_bottom_clearance.dispense = 16
        protocol_context.max_speeds['Z'] = 20 # Set Speed of Z Axis
        pipette.aspirate(aspirate_volume, dilution_buffer)
        protocol_context.max_speeds['Z'] = 200 # Set Speed of Z Axis
        pipette.dispense(aspirate_volume, target)
        protocol_context.delay(seconds=3)
        pipette.blow_out(target.bottom(16.5))
        pipette.touch_tip(v_offset=-5)
        protocol_context.max_speeds['Z'] = 200 # Set Speed of Z Axis
        pipette.drop_tip()
        Dilution_Buffer_Well_Bottom = Dilution_Buffer_Well_Bottom-4
        
    # Define Well Bottom for Reagent Reservoir (agilent_4_well_73_ml_reagent_reservoir) to dilution_10k_plate_plate (biorad_96_wellplate_200ul_pcr) 
    if col_num == 6:
        Dilution_Buffer_Well_Bottom=4
    if col_num == 5:
        Dilution_Buffer_Well_Bottom=3.4
    if col_num == 4:
        Dilution_Buffer_Well_Bottom=2.8
    if col_num == 3:
        Dilution_Buffer_Well_Bottom=2.2
    if col_num == 2:
        Dilution_Buffer_Well_Bottom=1.6
    if col_num == 1:
        Dilution_Buffer_Well_Bottom=1
     
    # Change Flow Rates
    pipette.flow_rate.aspirate = 35
    pipette.flow_rate.dispense = 95 # tested empircally, 95 gives better results than orignal set point of 250 JSB 08/30/21
    pipette.flow_rate.blow_out = 299
    dilution_volume_2 = (dilution_volume_2)-1 # giving a value of 95 to dispense actually yields 96 uL. so this is to correct the set volume to 94, which will actually yield 95 uL. 

    # Dispense Dilution Buffer to BioRad Hardshell 96-well plate (10k Dilution Plate)
    for target in samples_dilution_2:
        protocol_context.max_speeds['Z'] = 100 # Set Speed of Z Axis
        pipette.pick_up_tip()
        protocol_context.max_speeds['Z'] = 20 # Set Speed of Z Axis
        pipette.well_bottom_clearance.aspirate = Dilution_Buffer_Well_Bottom
        pipette.well_bottom_clearance.dispense = Dilution_Buffer_Well_Bottom
        pipette.mix(1, dilution_volume_2, dilution_buffer)
        pipette.well_bottom_clearance.aspirate = Dilution_Buffer_Well_Bottom
        pipette.well_bottom_clearance.dispense = 10
        pipette.aspirate(dilution_volume_2, dilution_buffer)
        protocol_context.max_speeds['Z'] = 100 # Set Speed of Z Axis
        pipette.dispense(dilution_volume_2, target)
        protocol_context.delay(seconds=3)
        pipette.blow_out(target.bottom(9.5))
        pipette.touch_tip(v_offset=-5)
        protocol_context.max_speeds['Z'] = 100 # Set Speed of Z Axis
        pipette.drop_tip()
        Dilution_Buffer_Well_Bottom = Dilution_Buffer_Well_Bottom-0.6
        
   # Dispense Dilution Buffer to 20k Dilution BioRad Hard Shell 96-well plate onto Temperature Module
   
   # Change Flow Rates
    pipette.flow_rate.aspirate = 35
    pipette.flow_rate.dispense = 95
    pipette.flow_rate.blow_out = 299
    dilution_volume_3 = (dilution_volume_3)-4 # giving a value of 36 to dispense actually yields 40-41 uL. so this is to correct the set volume to 39, which will actually yield 40-41 uL.
    
    Dilution_Buffer_Well_Bottom=1
    
    for target in samples_dilution_3:
        protocol_context.max_speeds['Z'] = 100 # Set Speed of Z Axis
        pipette.pick_up_tip()
        protocol_context.max_speeds['Z'] = 20 # Set Speed of Z Axis
        pipette.well_bottom_clearance.aspirate = Dilution_Buffer_Well_Bottom
        pipette.well_bottom_clearance.dispense = Dilution_Buffer_Well_Bottom
        pipette.mix(1, dilution_volume_3, dilution_buffer)
        pipette.well_bottom_clearance.aspirate = Dilution_Buffer_Well_Bottom
        pipette.well_bottom_clearance.dispense = 8
        pipette.aspirate(dilution_volume_3, dilution_buffer)
        protocol_context.max_speeds['Z'] = 100 # Set Speed of Z Axis
        pipette.dispense(dilution_volume_3, target)
        protocol_context.delay(seconds=3)
        pipette.blow_out(target.bottom(4))
        pipette.touch_tip(v_offset=-5)
        protocol_context.max_speeds['Z'] = 100 # Set Speed of Z Axis
        pipette.drop_tip()
        
    # Change Flow Rates
    pipette_2.flow_rate.aspirate = 2
    pipette_2.flow_rate.dispense = 2
    pipette_2.flow_rate.blow_out = 20
    #sample_volume_1 = (sample_volume_1)-0.2 # giving a value of 39 to dispense actually yields 40-41 uL. so this is to correct the set volume to 39, which will actually yield 40-41 uL.    
  
    #Initial Sample Transfer to PE Pipetting Microplate 2mL DW SQ 96-well plate
    for t in range(col_num):
        if t == 0:
            samples_1 = sample_plate['A1']
            dilutions_1 = dilution_500_plate['A1']
        if t == 1:
            samples_1 = sample_plate['A2']
            dilutions_1 = dilution_500_plate['A2']
        if t == 2:
            samples_1 = sample_plate['A3']
            dilutions_1 = dilution_500_plate['A3']
        if t == 3:
            samples_1 = sample_plate['A4']
            dilutions_1 = dilution_500_plate['A4']
        if t == 4:
            samples_1 = sample_plate['A5']
            dilutions_1 = dilution_500_plate['A5']
        if t == 5:
            samples_1 = sample_plate['A6']
            dilutions_1 = dilution_500_plate['A6']

        pipette_2.pick_up_tip()
        pipette_2.well_bottom_clearance.aspirate = starting_aspirate_height
        pipette_2.well_bottom_clearance.dispense = 16.5
        pipette_2.aspirate(sample_volume_1, samples_1)
        pipette_2.dispense(sample_volume_1, dilutions_1)
        pipette_2.well_bottom_clearance.aspirate = 16.5
        pipette_2.mix(1, sample_volume_1, dilutions_1)
        protocol_context.delay(seconds=3)
        pipette_2.blow_out(dilutions_1.bottom(17))
        pipette_2.touch_tip(v_offset=-5)
        pipette_2.drop_tip()
        
    protocol_context.set_rail_lights(True)
    
    #Set Temperature of Temperature Module to User Defined Variable
    temp_deck.set_temperature(set_temperature)
    
    #################################################################################################################################################
    ##########Prompt to Vortex PE Pipetting Microplate 2mL DW SQ 96-well plate############################################################################################################################################ 
    protocol_context.pause("Remove PE Pipetting Microplate 2mL DW SQ 96-well plate from deck grid 2. Seal Plate and Vortex(1 minute @ top speed)")
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("Centrifuge (*use DW bucket!) PE Pipetting Microplate 2mL DW SQ 96-well plate briefly until no bubbles remain (@200 x g, 1 minute). ReLoad PE Pipetting Microplate 2mL DW SQ 96-well plate onto deck grid 2. Remove seal. Once you click resume, pipetting will begin!")
    ##############################################################################################################################################################################
    ##############################################################################################################################################################################
   
    protocol_context.set_rail_lights(False)
    
    # Change Flow Rates
    pipette_2.flow_rate.aspirate = 5
    pipette_2.flow_rate.dispense = 5
    pipette_2.flow_rate.blow_out = 20

    # PE Pipetting Microplate 2mL DW SQ 96-well plate Sample Transfer to Dilution 10K Plate
    for t2 in range(col_num):
        if t2 == 0:
            samples_2 = dilution_500_plate['A1']
            dilutions_2 = dilution_10k_plate['A1']
        if t2 == 1:
            samples_2 = dilution_500_plate['A2']
            dilutions_2 = dilution_10k_plate['A2']
        if t2 == 2:
            samples_2 = dilution_500_plate['A3']
            dilutions_2 = dilution_10k_plate['A3']
        if t2 == 3:
            samples_2 = dilution_500_plate['A4']
            dilutions_2 = dilution_10k_plate['A4']
        if t2 == 4:
            samples_2 = dilution_500_plate['A5']
            dilutions_2 = dilution_10k_plate['A5']
        if t2 == 5:
            samples_2 = dilution_500_plate['A6']
            dilutions_2 = dilution_10k_plate['A6']

        pipette_2.pick_up_tip()
        pipette_2.well_bottom_clearance.aspirate = 15
        pipette_2.well_bottom_clearance.dispense = 8
        pipette_2.aspirate(sample_volume_2, samples_2)
        pipette_2.dispense(sample_volume_2, dilutions_2)
        pipette_2.well_bottom_clearance.aspirate = 8
        pipette_2.mix(1, sample_volume_2, dilutions_2)
        protocol_context.delay(seconds=3)
        pipette_2.blow_out(dilutions_2.bottom(9))
        pipette_2.touch_tip(v_offset=-5)
        pipette_2.drop_tip()
 
    protocol_context.set_rail_lights(True)
    
    #############################################################################################################################################################################
    ##########################################################Prompt to Vortex 10k Dilution BioRad Hard Shell 96-well plate###################################################
    #############################################################################################################################################################################  
    protocol_context.pause("Remove 10k Dilution BioRad Hard Shell 96-well plate from deck grid 3. Seal Plate and Vortex(1 minute @ top speed)")
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("Centrifuge plate briefly (@1500 rpm, 2 minutes). ReLoad 10k Dilution BioRad Hard Shell 96-well plate onto deck grid 3. Remove seal.Once you click resume, pipetting will begin!")
    ##############################################################################################################################################################################
    ##############################################################################################################################################################################
   
    protocol_context.set_rail_lights(False)

    # Change Flow Rates
    pipette.flow_rate.aspirate = 40
    pipette.flow_rate.dispense = 22.5
    pipette.flow_rate.blow_out = 299
    
    protocol_context.max_speeds['Z'] = 20 # Set Speed of(Z) Axis 20uL Pippete
          
    # Dilution 10K Plate Sample Transfer to Dilution 20K Plate
    for t3 in range(col_num):
        if t3 == 0:
            samples_3 = dilution_10k_plate['A1']
            dilutions_3 = dilution_20k_plate['A1']
        if t3 == 1:
            samples_3 = dilution_10k_plate['A2']
            dilutions_3 = dilution_20k_plate['A2']
        if t3 == 2:
            samples_3 = dilution_10k_plate['A3']
            dilutions_3 = dilution_20k_plate['A3']
        if t3 == 3:
            samples_3 = dilution_10k_plate['A4']
            dilutions_3 = dilution_20k_plate['A4']
        if t3 == 4:
            samples_3 = dilution_10k_plate['A5']
            dilutions_3 = dilution_20k_plate['A5']
        if t3 == 5:
            samples_3 = dilution_10k_plate['A6']
            dilutions_3 = dilution_20k_plate['A6']
        
        pipette.pick_up_tip()
        pipette.well_bottom_clearance.aspirate = 5
        pipette.well_bottom_clearance.dispense = 5
        pipette.aspirate(sample_volume_3, samples_3)
        pipette.dispense(sample_volume_3, dilutions_3)
        pipette.well_bottom_clearance.aspirate = 5
        pipette.mix(1, sample_volume_3, dilutions_3)
        protocol_context.delay(seconds=3)
        pipette.blow_out(dilutions_3.bottom(5.5))
        pipette.touch_tip(v_offset=-5)
        pipette.drop_tip()
        
    #Turn Lights On
    protocol_context.set_rail_lights(True)
    
    #############################################################################################################################################################################
    ##########################################################Prompt to Vortex 20k Dilution BioRad Hard Shell 96-well plate######################################################
    #############################################################################################################################################################################  
    protocol_context.pause("Remove 20k Dilution BioRad Hard Shell 96-well plate from Temperature Module on deck grid 10. Seal Plate and Vortex(1 minute @ top speed). Centrifuge Briefly(@1500 rpm, 2 minutes). Store on ice until ready to load onto deck for part 2 of the qPCR assay.")
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("Remove 10k Dilution BioRad Hard Shell 96-well plate from deck grid 3. Seal Plate and Vortex(1 minute @ top speed). Centrifuge Briefly(@1500 rpm, 2 minutes).Store on ice until ready to load onto deck for part 2 of the qPCR assay.")
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    ########################################################################Deck Unloading Instructions##########################################################################
    protocol_context.pause("Please remove Tip Waste from deck grid 12 to biohazard bin.")
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("Remove Initial Sample Plate from deck grid 1.Seal Plate, and store at 4C until qPCR data analysis is complete.")
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("Remove PE Pipetting Microplate 2mL DW SQ 96-well plate from deck grid 2 to biohazard bin.")
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("Remove Agilent 73 mL Reagent Reservoir to biohazard bin.")
    protocol_context.set_rail_lights(False)
    protocol_context.delay(seconds=1)
    protocol_context.set_rail_lights(True)
    protocol_context.pause("Remove Tip Boxes from deck. Wipe Down Deck with 70% Ethanol. Once you click 'Resume', pipette head will raise and script will end!")
    ##############################################################################################################################################################################
    ##############################################################################################################################################################################
    
     #Turn Lights Off, all Protocols Expect the lights to be off before starting
    protocol_context.set_rail_lights(False)
    
    # Turn off Temperature Deck, currently deactivated because you usually proceed to Step 2 Immediately
    #temp_deck.deactivate()
    
if __name__ == '__main__':
    protocol = simulate.get_protocol_api('2.10')
    run(protocol)
    for line in protocol.commands():
        print(line)
