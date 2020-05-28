"""
File: brickbreaker_protein_coding_edition.py
----------------
This program allows a user to play the game brick breaker.
It was written by Dia Bagchi in May 2020 as her Final Project
for Code in Place 2020.
"""

import tkinter
import time
import random

# How big is the playing area?
CANVAS_WIDTH = 900      # Width of drawing canvas in pixels
CANVAS_HEIGHT = 700     # Height of drawing canvas in pixels


# Constants for the bricks
N_ROWS = 8              # How many rows of bricks are there?
N_COLS = 15             # How many columns of bricks are there?
SPACING = 5             # How much space is there between each brick?
BRICK_START_Y = 50      # The y coordinate of the top-most brick
BRICK_HEIGHT = 20       # How many pixels high is each brick
BRICK_WIDTH = (CANVAS_WIDTH - (N_COLS+1) * SPACING ) / N_COLS
BRICK_COLORS = ["blue", "red", "green", "yellow2"]  # DNA colors
#BRICK_COLORS = ["burlywood3", "saddle brown", "cornsilk2", "OrangeRed3"]  #real brick colors
#BRICK_COLORS = ["blue", "magenta2", "lawn green", "ivory4"]  #fun colors

# Constants for the ball and paddle
BALL_SIZE = 40
BALL_COLOR = 'black'
PADDLE_Y = CANVAS_HEIGHT - 40
PADDLE_WIDTH = 80
PADDLE_HEIGHT = 10
PADDLE_COLOR = 'black'

# Constants for the ball speed
BALL_SPEED = 10
SLEEP_TIME = 1/50

# Constants for text
FONT = 'Helvetica'
FONT_SIZE = '30'
FONT_SIZE_2 = '20'

# Constants for protein blocks
PROTEIN_BLOCK_SIZE = 27
DNA_BLOCK_SIZE = 9
OFFSET = 2
BLOCK_START_Y = ((CANVAS_HEIGHT / 2) + 2 * (CANVAS_HEIGHT / 10))


def main():

    """
    Create a canvas to contain the game
    """
    canvas = make_canvas(CANVAS_WIDTH, CANVAS_HEIGHT, 'Brick Breaker')


    """
    Set up the world
    """
    # Make bricks
    make_bricks(canvas)

    # Introduce the rules of the game to the player
    introduction(canvas)

    # Make ball
    ball = make_ball(canvas)
    change_x = BALL_SPEED  # constant to define ball movement in x direction
    change_y = BALL_SPEED  # constant to define ball movement in y direction

    # Make paddle
    paddle = make_paddle(canvas)

    # Pause after initial set up to allow player to get ready
    move_just_paddle(canvas, paddle, 2)


    """
    Play the game
    """
    # Create a variable to count the number of turns that have gone by.
    game_turns = 0

    # Create a list to store the flavor of brick or dNTP hit
    nuc_list = []

    # Play the game while there are still bricks on the canvas
    while det_bricks_left(canvas, ball, paddle):
        # Count the turns
        game_turns = count_turns(canvas, ball, paddle, game_turns)
        # If the player has used all 3 turns, exit the game
        if game_turns >= 3:
            break

        """
        Move and monitor the ball
        """
        # Get the coordinates of the ball
        ball_x_coord = get_left_x(canvas, ball)
        ball_y_coord = get_top_y(canvas, ball)
        # If the ball hits the top wall change y direction to positive (down)
        if not (ball_y_coord > 0):
            change_y = BALL_SPEED
        # If the ball hits the paddle change y direction to negative (up)
        if paddle in canvas.find_overlapping(ball_x_coord + BALL_SIZE/4, ball_y_coord + BALL_SIZE/2, ball_x_coord + (BALL_SIZE*.75), ball_y_coord + BALL_SIZE):
            change_y = -BALL_SPEED
        # If the ball hits a side wall change x direction to opposite direction
        if not ( (ball_x_coord > 0) and (ball_x_coord < (CANVAS_WIDTH - BALL_SIZE)) ):
            change_x *= (-1)
        # Detect if the ball has hit any bricks
        # If it has, delete the bricks and change y direction to positive (down)
        bricks_hit = detect_hit_bricks(canvas, ball, change_y, paddle)
        if bricks_hit:
            nuc_list = flavor(canvas, bricks_hit, nuc_list)
            delete_bricks(canvas, bricks_hit)
            # If at least 1 brick was hit, the y direction of the ball should be positive
            change_y = BALL_SPEED
        canvas.move(ball, change_x, change_y)

        """
        Move and monitor the paddle
        """
        # Move Paddle to where the user places the mouse
        move_paddle(canvas, paddle)

        """
        Update the canvas after each loop and sleep for the set amount of time
        """
        canvas.update()
        time.sleep(SLEEP_TIME)


    """
    Evaluate the player's performance
    """
    conclusion(canvas, ball, paddle, nuc_list)


    canvas.mainloop()



def make_bricks(canvas):
    """
    This function produces all the bricks.
    """
    brick_list = []
    for row in range(N_ROWS):
        for col in range(N_COLS):
            brick_x1 = ((SPACING + BRICK_WIDTH) * col) + SPACING
            brick_x2 = (SPACING + BRICK_WIDTH) * (col + 1)
            brick_y1 = BRICK_START_Y + ((SPACING + BRICK_HEIGHT) * row)
            brick_y2 = BRICK_START_Y + (BRICK_HEIGHT * (row + 1)) + (SPACING * row)
            color = brick_colors()
            brick_list.append(canvas.create_rectangle(brick_x1, brick_y1, brick_x2, brick_y2, fill=color))

    return brick_list


def brick_colors():
    """
    This function pulls a random color out of a list of colors.
    The color list is a global variable.
    """
    colors = BRICK_COLORS

    return colors[random.randint(0,3)]


def make_ball(canvas):
    """
    This function creates the ball object on the canvas.
    """
    ball_x1 = (CANVAS_WIDTH / 2) - (BALL_SIZE / 2)
    ball_x2 = ball_x1 + BALL_SIZE
    ball_y1 = (CANVAS_HEIGHT / 2) - (BALL_SIZE / 2)
    ball_y2 = ball_y1 + BALL_SIZE
    ball = canvas.create_oval(ball_x1, ball_y1, ball_x2, ball_y2, fill=BALL_COLOR)

    return ball


def make_paddle(canvas):
    """
    This function creates the paddle object on the canvas.
    """
    paddle = canvas.create_rectangle(5, PADDLE_Y, (PADDLE_WIDTH +5), (PADDLE_Y + PADDLE_HEIGHT), fill=PADDLE_COLOR)

    return paddle


def introduction(canvas):
    """
    This function displays the game title page and gives
    the player instructions for how to play the game.
    """

    # Produce the game's colorful title page
    # Display the game title (colored letters)
    text_list = []  # Variable to store text objects so that they may be deleted afterwards
    message = "Hello!"
    text_list.extend(print_fun_text(canvas, message, 1))
    message = "Welcome to Brickbreaker!"
    text_list.extend(print_fun_text(canvas, message, 2))
    time.sleep(1)

    # Display the game subtitle (colored words)
    message = "Protein"
    text_list.append(canvas.create_text((CANVAS_WIDTH*.4), ((CANVAS_HEIGHT / 2) + 2.5 * (CANVAS_HEIGHT / 10)),
                                        text=message, anchor='center', font=(FONT, FONT_SIZE_2), fill=brick_colors()))
    canvas.update()
    time.sleep(.75)
    message = "Coding"
    text_list.append(canvas.create_text((CANVAS_WIDTH / 2), ((CANVAS_HEIGHT / 2) + 2.5 * (CANVAS_HEIGHT / 10)),
                                        text=message, anchor='center', font=(FONT, FONT_SIZE_2), fill=brick_colors()))
    canvas.update()
    time.sleep(.75)
    message = "Edition"
    text_list.append(canvas.create_text((CANVAS_WIDTH*.6), ((CANVAS_HEIGHT / 2) + 2.5 * (CANVAS_HEIGHT / 10)),
                                        text=message, anchor='center', font=(FONT, FONT_SIZE_2), fill=brick_colors()))
    canvas.update()
    time.sleep(1)

    # Display attribution info
    message = "Created by Dia Bagchi"
    text_list.append(canvas.create_text((CANVAS_WIDTH / 2), ((CANVAS_HEIGHT / 2) + 4 * (CANVAS_HEIGHT/10)),
                                      text=message, anchor='center', font=(FONT, FONT_SIZE_2)))
    canvas.update()
    time.sleep(1.5)
    message = "for Code in Place 2020 Final Project"
    text_list.append(canvas.create_text((CANVAS_WIDTH / 2), ((CANVAS_HEIGHT / 2) + 4.5 * (CANVAS_HEIGHT / 10)),
                                        text=message, anchor='center', font=(FONT, FONT_SIZE_2)))
    canvas.update()
    time.sleep(4)

    # Delete all the text from the canvas
    delete_objects(canvas, text_list)


    # Give player instructions
    # Define positions for multiple lines of text
    line_1_y = ((CANVAS_HEIGHT / 2) + (CANVAS_HEIGHT / 10))
    line_2_y = ((CANVAS_HEIGHT / 2) + (CANVAS_HEIGHT / 10)*1.5)
    line_3_y = ((CANVAS_HEIGHT / 2) + (CANVAS_HEIGHT / 10) * 2)
    line_4_y = ((CANVAS_HEIGHT / 2) + (CANVAS_HEIGHT / 10) * 2.5)

    # Display instructions
    message = "Use the paddle to keep the ball bouncing and break all the bricks."
    display_instructions(canvas, message, line_1_y, secs=4)

    message = "You will have 3 chances to hit all the bricks before the game ends."
    display_instructions(canvas, message, line_1_y, secs=4)

    message = "In this version of the game"
    display_instructions(canvas, message, line_1_y, secs=2.5)

    message = "the bricks have been color coded to represent DNA nucleotides."
    display_instructions(canvas, message, line_1_y, secs=4)


    # Display the nucleotide color coding scheme
    A_nuc = "Adenine is red"
    T_nuc = "Thymine is yellow"
    G_nuc = "Guanine is blue"
    C_nuc = "Cytosine is green"
    text_list = []

    text_list.append(display_nuc_coding(canvas, A_nuc, line_1_y, 'red'))
    text_list.append(display_nuc_coding(canvas, T_nuc, line_2_y, 'yellow2'))
    text_list.append(display_nuc_coding(canvas, G_nuc, line_3_y, 'blue'))
    text_list.append(display_nuc_coding(canvas, C_nuc, line_4_y, 'green'))
    time.sleep(3)
    # Delete the color coding scheme
    delete_objects(canvas, text_list)


    # Continue instructions
    message = "DNA sequences which encode a protein can be found when"
    display_instructions(canvas, message, line_1_y, secs=4.5)


    # Explanation of start codon
    text_list = []
    # Left side text
    message = "the 3 nucleotides"
    text_list.append(canvas.create_text((CANVAS_WIDTH * .3), line_1_y,
                                        text=message, anchor='center', font=(FONT, FONT_SIZE_2) ))
    # Center text -display colored letters for the nucleotides found in the start codon
    text_list.append(canvas.create_text((CANVAS_WIDTH * .45), line_1_y,
                                        text="A", anchor='center', font=(FONT, FONT_SIZE), fill='red'))
    text_list.append(canvas.create_text((CANVAS_WIDTH * .5), line_1_y,
                                        text="T", anchor='center', font=(FONT, FONT_SIZE), fill='yellow2'))
    text_list.append(canvas.create_text((CANVAS_WIDTH * .55), line_1_y,
                                        text="G", anchor='center', font=(FONT, FONT_SIZE), fill='blue'))
    # Right side text
    message = "appear sequentially."
    text_list.append(canvas.create_text((CANVAS_WIDTH * .7), line_1_y,
                                        text=message, anchor='center', font=(FONT, FONT_SIZE_2), ))
    canvas.update()
    time.sleep(5)
    delete_objects(canvas, text_list)


    # Continue instructions
    message = "This sequence is called a start codon, and can mark the start of a protein coding sequence in DNA."
    display_instructions(canvas, message, line_1_y, secs=6)

    message = "If the sequence in which you break the bricks encodes possible proteins,"
    display_instructions(canvas, message, line_1_y, secs=5)

    message = "you will be shown the translated sequences at the end of the game."
    display_instructions(canvas, message, line_1_y, secs=5)



def print_fun_text(canvas, message, text_line):
    """
    This function displays text in a colorful manner.  Each letter
    in the words in the input text is assigned a random color from
    a list of allowed colors.  The letters are then alternately offset
    from each other.
    """
    if text_line == 1:
        text_start_y_pos = (CANVAS_HEIGHT / 2)
    if text_line == 2:
        text_start_y_pos = (CANVAS_HEIGHT / 2) + 1.25*(CANVAS_HEIGHT /10)
    else:
        text_start_y_pos = (CANVAS_HEIGHT / 2)

    spacing = 20
    text_width = ((len(message) *spacing))
    text_start = (CANVAS_WIDTH/2)-(text_width/2)
    text_list = []
    count = 0
    for i in range(len(message)):
        if count % 2 == 0:
            text_y_pos = text_start_y_pos -3
        else:
            text_y_pos = text_start_y_pos +3
        if message[i] != " ":
            letter = canvas.create_text((text_start + spacing * i), text_y_pos, text=message[i], anchor='center', font=(FONT, FONT_SIZE), fill=brick_colors())
            text_list.append(letter)
            count += 1
        canvas.update()
        # pause
        time.sleep(.2)

    return text_list


def display_instructions(canvas, message, height, secs):
    """
    This function formats the display of player instructions on the canvas.
    """
    screen_text = canvas.create_text((CANVAS_WIDTH / 2), height,
                                     text=message, anchor='center', font=(FONT, FONT_SIZE_2), fill='black')
    canvas.update()
    # pause
    time.sleep(secs)
    canvas.delete(screen_text)


def display_nuc_coding(canvas, message, height, color):
    """
    This function displays colored text to demonstrate
    the nucleotide color coding scheme.
    """
    screen_text = canvas.create_text((CANVAS_WIDTH / 2), height,
                                     text=message, anchor='center', font=(FONT, FONT_SIZE_2), fill=color)
    canvas.update()
    # pause
    time.sleep(1.75)
    return screen_text


def pause(canvas, secs):
    """
    This function updates the canvas and pauses for a set amount of time.
    """
    canvas.update()
    # pause
    time.sleep(secs)


def move_paddle(canvas, paddle):
    """
    This function controls the movement of the paddle.
    The paddle moves along the x axis based on the location
    of the cursor along the axis.  The paddle has been adjusted
    to be centered on the cursor position.  The y axis location
    of the paddle is fixed, meaning that the paddle can only
    move from left to right.
    """
    mouse_x = canvas.winfo_pointerx()
    paddle_x_centered = mouse_x - (PADDLE_WIDTH/2)
    if paddle_x_centered > CANVAS_WIDTH - (PADDLE_WIDTH/2):
        paddle_x_centered = CANVAS_WIDTH - PADDLE_WIDTH
    if paddle_x_centered < 0 + (PADDLE_WIDTH/2):
        paddle_x_centered = 0
    canvas.moveto(paddle, paddle_x_centered, PADDLE_Y)


def det_bricks_left(canvas, ball, paddle):
    """
    This function determines how many bricks are left.
    """
    bricks = []
    canvas_objects = canvas.find_overlapping(0, 0, CANVAS_WIDTH, CANVAS_HEIGHT)
    for elem in canvas_objects:
        # Test if the element is not the ball or paddle
        if elem != paddle and elem != ball:  # Assumes all other elements are bricks
            bricks.append(elem)
    return bricks


def ball_overlaps(canvas, ball):
    """
    This function finds the coordinates on the Canvas taken up by the ball.
    It then checks for overlapping objects, adds any found to a list, and
    returns the list.
    """
    ball_coords = canvas.coords(ball)  # returns a list with 4 coordinates
    x_1 = ball_coords[0]
    y_1 = ball_coords[1]
    x_2 = ball_coords[2]
    y_2 = ball_coords[3]
    colliding_list = canvas.find_overlapping(x_1, y_1, x_2, y_2)
    return colliding_list


def detect_hit_bricks(canvas, ball, change_y, paddle):
    """
    This function determines if any bricks were hit by the ball.
    It then returns a list of the bricks hit.
    """
    # Find ball overlaps
    colliding_list = ball_overlaps(canvas, ball)
    # Start by assuming no bricks overlap the ball
    bricks_hit = []
    brick = False
    # Test whether the ball hit a brick
    if len(colliding_list) >= 2:  # Find if anything overlaps the ball
        # For each overlapping element
        for elem in colliding_list:
            # Test if the element is not the ball or paddle
            if elem != paddle and elem != ball:  # Assumes all other elements are bricks
                bricks_hit.append(elem)
    return bricks_hit


def flavor(canvas, bricks_hit, nuc_list):
    """
    This function converts the color of any brick hit
    into the corresponding DNA nucleotide.
    """
    for item in bricks_hit:
        flavour = canvas.itemcget(item, 'fill')
        if flavour == "red":
            nucleotide = "A"
        if flavour == "yellow2":
            nucleotide = "T"
        if flavour == "blue":
            nucleotide = "G"
        if flavour == "green":
            nucleotide = "C"
    nuc_list.append(nucleotide)

    return nuc_list


def delete_bricks(canvas, bricks):
    """
    This function deletes hit bricks.
    """
    # For each brick in the list
    for brick in bricks:
        # Delete the brick
        canvas.delete(brick)


def count_turns(canvas, ball, paddle, game_turns):
    """
    This function counts the number of turns that the player has had.
    """
    if ball not in canvas.find_overlapping(0, 0, CANVAS_WIDTH, CANVAS_HEIGHT):
        # Move ball completely out of the canvas
        canvas.moveto(ball, (CANVAS_WIDTH / 2), (CANVAS_HEIGHT * 2))
        game_turns += 1
        if game_turns < 3:
            print_turns(canvas, paddle, game_turns)
            canvas.moveto(ball, (CANVAS_WIDTH / 2 - BALL_SIZE / 2), (CANVAS_HEIGHT / 2 - BALL_SIZE / 2))

            move_just_paddle(canvas, paddle, 1)

    return game_turns


def print_turns(canvas, paddle, game_turns):
    """
    This function displays the number of turns a player has left
    once they have missed hitting the ball.
    """
    message = "Try again: " + str(3 - game_turns) + " out of 3 chances left."
    game_text = canvas.create_text((CANVAS_WIDTH / 2), (CANVAS_HEIGHT / 2),
                                   text=message, anchor='center', font=(FONT, FONT_SIZE))
    move_just_paddle(canvas, paddle, 2)
    canvas.delete(game_text)


def move_just_paddle(canvas, paddle, secs):
    """
    This function allows the paddle to be moved while the
    rest of the canvas is held static for a set amount of time.
    This allows the player to get ready for the next turn.
    """
    start_time = time.time()
    while time.time() < start_time + secs:
        move_paddle(canvas, paddle)
        canvas.update()


def conclusion(canvas, ball, paddle, nuc_list):
    """
    This function evaluates how the player performed, and
    reports back if any open reading frames (ORFs) have been
    discovered in the sequence of bricks broken.  It also
    translates the ORF DNA sequence into the corresponding
    amino acids they code for.
    """
    # Evaluate player performance with regard to breaking all the bricks.
    if det_bricks_left(canvas, ball, paddle):
        message = "Game over."
        canvas.create_text((CANVAS_WIDTH / 2), (CANVAS_HEIGHT / 2), text=message, anchor='center',
                           font=(FONT, FONT_SIZE))
    else:
        message = "Congratulations!"
        message2 = "There are no more bricks!"
        canvas.create_text((CANVAS_WIDTH / 2), (CANVAS_HEIGHT / 2), text=message, anchor='center',
                           font=(FONT, FONT_SIZE))
        canvas.create_text((CANVAS_WIDTH / 2), ((CANVAS_HEIGHT / 2) + (CANVAS_HEIGHT/10)), text=message2, anchor='center',
                           font=(FONT, FONT_SIZE_2))
    # Search for a potential ORF in the order of bricks broken.
    orfs_list = find_orfs(nuc_list)
    # If there is an ORF
    if orfs_list:
        # Translate the nucleotide sequence into a protein sequence, and display the sequences as blocks.
        if len(orfs_list)> 3:
            orfs_list = truncate_orfs_list(orfs_list)
        count = 0
        for orf in orfs_list:
            produce_protein(canvas, orf, count)
            count += 1
        if len(orfs_list) == 1:
            message = "Possible protein sequence found!"
        else:
            message = "Possible protein sequences found!"
        canvas.create_text((CANVAS_WIDTH / 2), ((CANVAS_HEIGHT / 2) + 4 * (CANVAS_HEIGHT / 10)), text=message,
                           anchor='center',
                           font=(FONT, FONT_SIZE_2))

    # If no ORF was found
    else:
        message = "No potential coding sequences were found in the series of bricks broken."
        canvas.create_text((CANVAS_WIDTH / 2), ((CANVAS_HEIGHT / 2) + 4 * (CANVAS_HEIGHT / 10)), text=message,
                           anchor='center',
                           font=(FONT, FONT_SIZE_2))
    canvas.update()
    # Allow paddle to move for 10 secs after the game is over.
    move_just_paddle(canvas, paddle, 10)
    # Delete the paddle.
    canvas.delete(paddle)



def find_orfs(nuc_list):
    """
    This function takes in a list of nucleotides that correspond
    to the order (by color) that the bricks were broken in. It then
    searches the list for a start codon (ATG).  If a start codon is
    found, then the codons of the open reading frame are gathered until
    a stop codon is found or the sequence ends.
    """
    orfs_list = []  # Variable to store the ORF sequence
    start_index = None  # Assume there are no start codons to begin with
    start_indices_list = []
    for i in range(len(nuc_list)-2):
        triplet_nuc = nuc_list[i:i+3]  # Read the nucleotides as codons (groups of 3)
        if triplet_nuc == ['A', 'T', 'G']:  # If a start codon is found
            start_index = i  # Save the index of its first nucleotide (A)
            start_indices_list.append(start_index)
    if start_indices_list:  # If a start codon was identified
        for index in start_indices_list:

            coding_orf_nucs_list = []
            for i in range(index, len(nuc_list)-2, 3):
                codon = nuc_list[i:i+3]  # Grab codons starting from its reading frame
                if len(codon) == 3:
                    if codon != ['T', 'A', 'A'] and codon != ['T', 'G', 'A'] and codon != ['T', 'A', 'G']:
                        coding_orf_nucs_list.extend(codon)  # Keep reading codons if there is no stop
                    else:  # If there is a stop codon
                        coding_orf_nucs_list.extend(codon)  # Add the stop codon to the sequence
                        break  # Break out of the loop (do not add more nucleotides after stop)
            orfs_list.append(coding_orf_nucs_list)

    return orfs_list  # Return coding nucleotides


def truncate_orfs_list(orfs_list):
    """
    This function limits the number of potential protein sequences
    displayed at the end of the game to 3.  If more than 3 ORFs were
    detected, only the 3 longest ones are displayed for the player.
    """
    orfs_list.sort(key=len)
    trunc_list = [orfs_list[-1], orfs_list[-2], orfs_list[-3]]
    return trunc_list


def produce_protein(canvas, orf_list, count):
    """
    This function converts the identified coding nucleotide seq
    to a protein sequence.  It then produces blocks to represent them.
    """
    dna_seq = ""
    for elem in orf_list:
        dna_seq += str(elem)
    # Convert the coding sequence to an amino sequence
    protein = make_protein_seq(dna_seq)
    # Produce the blocks to represent the protein amino acids and DNA nucleotides
    create_blocks(canvas, protein, dna_seq, count)


def make_protein_seq(dna_seq):
    """
    This function takes in a DNA sequence and returns the corresponding
    protein sequence.
    """
    protein = ""
    for i in range(0,len(dna_seq),3):
        codon = dna_seq[i:i + 3]
        amino_acid = codon_look_up(codon)
        protein += amino_acid
    return protein


def codon_look_up(codon):
    """
    This function takes in a codon, looks up the
    corresponding amino acid in a Codon Table, and returns the
    encoded amino acid. The Codon Table is in the form of a dictionary.
    """
    # This codon table was obtained from: https://www.geeksforgeeks.org/dna-protein-python-3/
    codon_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }

    return codon_table[codon]


def create_blocks(canvas, protein_seq, dna_seqs, count):
    """
    This function takes in protein and DNA sequences, and creates a block
    for each amino acid and nucleotide.
    """
    # The starting y position for the blocks is increased for each sequence displayed
    position_y = BLOCK_START_Y + count*(PROTEIN_BLOCK_SIZE + DNA_BLOCK_SIZE + OFFSET)
    position_x = (CANVAS_WIDTH / 2) - ((PROTEIN_BLOCK_SIZE + OFFSET) * (len(protein_seq) / 2))

    for amino_acid in protein_seq:
        color = color_look_up(amino_acid)
        create_block(canvas, amino_acid, position_x, position_y, color, "amino_acid")
        position_x += PROTEIN_BLOCK_SIZE + OFFSET
    position_y += PROTEIN_BLOCK_SIZE
    position_x = (CANVAS_WIDTH/2) - ((PROTEIN_BLOCK_SIZE + OFFSET) * (len(protein_seq)/2))
    nuc_count = 0
    for dntp in dna_seqs:
        if nuc_count % 3 == 0 and nuc_count != 0:
            position_x += OFFSET
        nuc_count += 1
        color = dntp_color_look_up(dntp)
        create_block(canvas, dntp, position_x, position_y, color, "dNTP")
        position_x += DNA_BLOCK_SIZE


def color_look_up(amino_acid):
    """
    This function takes in an amino acid, looks up the corresponding
    color in a table, and returns the associated color.
    """
    # Source for color scheme: https://pubmed.ncbi.nlm.nih.gov/22607364/
    color_table = {
    'P':'blue', 'N': 'dodger blue', 'D':'deep sky blue', 'S':'cyan',
    'Q':'forest green', 'E':'lime green', 'R':'lawn green',
    'M':'red', 'I':'orange red', 'L':'orange', 'F':'yellow',
    'Y':'tan4', 'W':'dark goldenrod', 'K':'goldenrod3', 'H':'goldenrod2',
    'A':'purple3', 'G':'DarkOrchid2', 'V':'magenta2', 'C':'deep pink',
    'T':'gray64', '_':'white'
    }

    color = color_table[amino_acid]
    return color


def dntp_color_look_up(dntp):
    """
    This function takes in a nucleotide, looks up the corresponding color
    in a table, and returns the associated color.
    """
    # Source for color scheme: https://pubmed.ncbi.nlm.nih.gov/22607364/
    color_table = {
    'A':'red', 'T':'yellow2', 'G': 'blue', 'C':'green',
    }

    color = color_table[dntp]
    return color


def create_block(canvas, text, position_x, position_y, color, type):
    """
    This function makes the individual blocks that get displayed.
    """
    if type == "amino_acid":
        block_size = PROTEIN_BLOCK_SIZE
        text_size = PROTEIN_BLOCK_SIZE
    else:
        block_size = DNA_BLOCK_SIZE
        text_size = DNA_BLOCK_SIZE
    canvas.create_rectangle(position_x, position_y, position_x + block_size, position_y + block_size, fill= color)
    canvas.create_text(position_x + block_size/2, position_y + block_size/2, text=str(text), anchor='center', font=('Courrier', text_size))


def delete_objects(canvas, obj_list):
    """
    This function deletes items in a list from the canvas.
    """
    for elem in obj_list:
        canvas.delete(elem)


def get_top_y(canvas, object):
    '''
    This friendly method returns the y coordinate of the top of an object.
    Recall that canvas.coords(object) returns a list of the object 
    bounding box: [x_1, y_1, x_2, y_2]. The element at index 1 is the top-y
    '''
    return canvas.coords(object)[1]

def get_left_x(canvas, object):
    '''
    This friendly method returns the x coordinate of the left of an object.
    Recall that canvas.coords(object) returns a list of the object 
    bounding box: [x_1, y_1, x_2, y_2]. The element at index 0 is the left-x
    '''
    return canvas.coords(object)[0]

def make_canvas(width, height, title):
    """
    DO NOT MODIFY
    Creates and returns a drawing canvas
    of the given int size with a blue border,
    ready for drawing.
    """
    top = tkinter.Tk()
    top.minsize(width=width, height=height)
    top.title(title)
    canvas = tkinter.Canvas(top, width=width + 1, height=height + 1)
    canvas.pack()
    return canvas

if __name__ == '__main__':
    main()
