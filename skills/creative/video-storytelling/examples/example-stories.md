# Example Stories

## Complete Story Walkthrough Examples

This document provides complete examples of video storytelling from planning through execution. Each example shows the full story structure, character definitions, scene descriptions, and generated outputs.

## Example 1: Pyter's First Bug

A lighthearted educational story about a Python snake programmer encountering and debugging his first bug.

### Story Metadata

```python
story = {
    "title": "Pyter's First Bug",
    "theme": "Problem-solving and friendship",
    "tone": "Lighthearted educational",
    "target_audience": "All ages",
    "duration": "~100 seconds (6 scenes)",
    "moral": "Sometimes the best debugging tool is a fresh perspective"
}
```

### Characters

```python
characters = {
    "Pyter": {
        "species": "python snake",
        "personality": "curious, friendly, determined",
        "goal": "write perfect code",
        "voice_attributes": {
            "role": "narrator",
            "age": "young",
            "gender": "male"
        },
        "voice_id": "JBFqnCBsd6RMkjVDRZzb",  # George - narrative voice
        "appearance": """
            Species: cartoon python snake
            Colors: bright blue body (#0B5FFF) with yellow accents (#FFB703)
            Features: round friendly eyes, small smile, coiled body showing 2-3 coils
            Clothing: small red developer beanie hat
            Size: medium-sized snake, about 3 feet visible
            Style: semi-realistic cartoon matching global style
        """
    },

    "Bug": {
        "species": "small beetle",
        "personality": "helpful, cheerful, innocent",
        "purpose": "teach Pyter about debugging",
        "appearance": """
            Species: small cartoon bug/beetle
            Colors: orange-red body (#FB8500) with dark spots (#023047)
            Features: large expressive eyes, tiny antennae, 6 legs
            Clothing: none
            Size: very small, fits on keyboard
            Style: semi-realistic cartoon matching global style
        """
    }
}
```

### Story Arc

```python
arc = [
    {"scene": 0, "beat": "Introduction", "tension": 0, "emotion": "neutral"},
    {"scene": 1, "beat": "Normal world", "tension": 1, "emotion": "happy"},
    {"scene": 2, "beat": "Inciting incident", "tension": 3, "emotion": "surprised"},
    {"scene": 3, "beat": "Investigation", "tension": 5, "emotion": "curious"},
    {"scene": 4, "beat": "Breakthrough", "tension": 7, "emotion": "excited"},
    {"scene": 5, "beat": "Resolution", "tension": 2, "emotion": "triumphant"}
]
```

### Complete Scenes

**Scene 0: Title - Meet Pyter**

```python
scene_0 = {
    "visual_description": """
        Pyter the Python snake sits at a modern computer desk, coiled comfortably
        around an ergonomic chair. His bright blue body contrasts beautifully with
        the warm desk lighting. He wears a small red developer beanie tilted slightly
        on his head. The computer monitor glows softly in front of him, displaying
        colorful code. His round friendly eyes look directly at the viewer with an
        welcoming expression. The background shows a clean, minimalist home office
        with subtle gradient from light blue to cream.
    """,

    "narration": """
        Meet Pyter, a Python snake who loves to code. Today he's working on his
        most ambitious project yet, a program to help other snakes learn programming.
        Little does he know, an unexpected visitor is about to change everything.
    """,

    "emotion": "neutral",
    "word_count": 58,
    "expected_duration": "15 seconds"
}
```

**Scene 1: Coding Time**

```python
scene_1 = {
    "visual_description": """
        Pyter is actively typing at his keyboard, his tail wrapped around the mouse.
        The monitor shows beautiful, well-formatted Python code scrolling smoothly.
        His expression is focused and happy, with a slight smile. The screen glows
        warmly, illuminating his face. Code syntax highlighting shows in vibrant
        colors matching the overall palette. Everything feels productive and peaceful.
        Same desk setup as before, maintaining visual consistency.
    """,

    "narration": """
        Pyter's code flowed smoothly, each line more elegant than the last. Functions
        nested perfectly, variables named just right. This was his happy place, where
        logic and creativity danced together in perfect harmony. Everything was going
        exactly according to plan.
    """,

    "emotion": "happy",
    "word_count": 72,
    "expected_duration": "18 seconds",
    "reference_previous": True
}
```

**Scene 2: The Bug Appears**

```python
scene_2 = {
    "visual_description": """
        Pyter's eyes suddenly widen in surprise. His head tilts down toward the
        keyboard. There, sitting right on the spacebar, is a small orange beetle
        with dark spots. The bug looks cheerful and innocent, unaware of Pyter's
        surprise. Pyter's expression shifts from focused to surprised, mouth slightly
        open. The lighting remains consistent, but focus shifts to highlight both
        Pyter's face and the small bug. Same desk and setup as previous scenes.
    """,

    "narration": """
        Pyter was deep in his code when something caught his eye. A tiny orange bug
        sat right on his keyboard! He'd seen bugs in his code before, but never a
        real bug on his computer. This was definitely not part of the documentation.
        What was going on?
    """,

    "emotion": "surprised",
    "word_count": 65,
    "expected_duration": "17 seconds",
    "reference_previous": True
}
```

**Scene 3: Investigation**

```python
scene_3 = {
    "visual_description": """
        Pyter leans in closer to examine the bug, his face now very near the keyboard.
        His expression is curious and focused, eyes slightly squinted in concentration.
        The orange bug sits calmly, looking back at Pyter with large, friendly eyes.
        Pyter's beanie tilts forward slightly from leaning. The computer screen still
        glows in the background, code visible but slightly out of focus. Close-up
        composition showing both Pyter's curious face and the bug clearly.
    """,

    "narration": """
        Pyter leaned in closer, studying the little creature. It had six legs and
        bright orange spots. Wait a minute... orange spots? Just like the error
        messages in his console! Could this bug be trying to help him find the bugs
        in his code? That would be remarkable!
    """,

    "emotion": "curious",
    "word_count": 68,
    "expected_duration": "17 seconds",
    "reference_previous": True
}
```

**Scene 4: The Fix**

```python
scene_4 = {
    "visual_description": """
        Pyter sits upright again, eyes bright with excitement. One of his coils
        points at the screen where line 42 of code is highlighted. The bug now sits
        on Pyter's tail, visible in the frame. Pyter's expression is excited and
        triumphant - he's found it! The screen shows a before/after comparison: a
        small typo circled in red, then the same line corrected in green. The bug
        looks equally happy, antennae perked up cheerfully.
    """,

    "narration": """
        Together they found it! Line forty-two had a tiny typo that Pyter had missed
        a dozen times. The bug's cheerful presence had helped him see it with fresh
        eyes. Sometimes the best debugging tool is a new perspective and a friend who
        cares about your success!
    """,

    "emotion": "excited",
    "word_count": 70,
    "expected_duration": "18 seconds",
    "reference_previous": True
}
```

**Scene 5: Success**

```python
scene_5 = {
    "visual_description": """
        The screen now displays "SUCCESS!" in large, bright green letters with
        cheerful checkmarks. Pyter's expression is pure joy - wide smile, bright
        eyes, victorious pose with one coil raised in celebration. The bug sits
        on Pyter's beanie, both of them celebrating together. The entire scene
        feels warm and triumphant. Code in the background shows all green status
        indicators. The lighting seems even warmer and more positive than before.
    """,

    "narration": """
        The program ran flawlessly! Pyter's code lit up the screen in beautiful
        success messages. He smiled at his new friend the bug, who had taught him
        that sometimes the best way to find bugs in your code is to find a real
        bug who cares. Mission accomplished!
    """,

    "emotion": "triumphant",
    "word_count": 67,
    "expected_duration": "17 seconds",
    "reference_previous": True
}
```

### Generation Commands

```python
# Full implementation
from pathlib import Path
import os
from openai import OpenAI
from elevenlabs.client import ElevenLabs

# Initialize clients
openai_client = OpenAI(api_key=os.environ["OPENAI_API_KEY"])
elevenlabs_client = ElevenLabs(api_key=os.environ["ELEVENLABS_API_KEY"])

# Style locks (as defined in skill)
STYLE_LOCK = """
Aspect ratio: 1080×1080 (square)
Camera: 50mm lens, eye-level perspective
Lighting: soft three-point lighting, warm key light (4500K)
Color palette: #0B5FFF, #FFB703, #FB8500, #023047, #8ECAE6
Materials: matte finish, no film grain or heavy bloom
Background: subtle gradient, clean composition
Style: semi-realistic cartoon with clear lines and gentle shading
Post: crisp focus, no vignette or text artifacts
"""

NEGATIVE_LOCK = """
Avoid: photorealism, heavy grain, excessive bloom, text overlays,
watermarks, vignette, lens distortion, blurry focus, cluttered
backgrounds, inconsistent lighting, mixed art styles
"""

# Generate all scenes
output_dir = Path("output/pyters_first_bug")
output_dir.mkdir(parents=True, exist_ok=True)

scenes = [scene_0, scene_1, scene_2, scene_3, scene_4, scene_5]
previous_image = None

for i, scene in enumerate(scenes):
    print(f"\nGenerating Scene {i}...")

    # Build prompt
    prompt = f"{STYLE_LOCK}\n\n{characters['Pyter']['appearance']}\n\n"
    if i >= 2:  # Bug appears in scene 2
        prompt += f"{characters['Bug']['appearance']}\n\n"

    if previous_image:
        prompt = f"REFERENCE IMAGE: Maintain exact visual style.\n\n{prompt}"

    prompt += f"Scene: {scene['visual_description']}\n\n{NEGATIVE_LOCK}"

    # Generate image
    if previous_image:
        # Use DALL-E 2 with variation for consistency
        response = openai_client.images.create_variation(
            image=open(previous_image, "rb"),
            n=1,
            size="1024x1024"
        )
    else:
        # First scene, use DALL-E 3
        response = openai_client.images.generate(
            model="dall-e-3",
            prompt=prompt,
            size="1024x1024",
            quality="standard",
            n=1
        )

    # Save image
    image_url = response.data[0].url
    image_path = output_dir / f"scene_{i:03d}.png"
    # Download and save (implementation detail)
    # save_image_from_url(image_url, image_path)
    previous_image = image_path

    # Generate audio
    audio = elevenlabs_client.text_to_speech.convert(
        text=scene['narration'],
        voice_id=characters['Pyter']['voice_id'],
        model_id="eleven_multilingual_v2",
        output_format="mp3_44100_128"
    )

    audio_path = output_dir / f"scene_{i:03d}.mp3"
    with open(audio_path, "wb") as f:
        for chunk in audio:
            f.write(chunk)

    print(f"  ✓ Image: {image_path}")
    print(f"  ✓ Audio: {audio_path}")

# Assemble video
print("\nAssembling video...")
subprocess.run([
    "bash",
    "scripts/assemble_video.sh",
    str(output_dir),
    "pyters_first_bug.mp4"
])

print(f"\n✅ Complete! Video: {output_dir}/pyters_first_bug.mp4")
```

### Expected Output

**Files Generated:**
```
output/pyters_first_bug/
├── scene_000.png          (Title: Pyter at desk)
├── scene_000.mp3          (15 seconds, neutral narration)
├── scene_001.png          (Pyter coding happily)
├── scene_001.mp3          (18 seconds, happy narration)
├── scene_002.png          (Bug appears!)
├── scene_002.mp3          (17 seconds, surprised narration)
├── scene_003.png          (Investigating the bug)
├── scene_003.mp3          (17 seconds, curious narration)
├── scene_004.png          (Finding the solution)
├── scene_004.mp3          (18 seconds, excited narration)
├── scene_005.png          (Success celebration)
├── scene_005.mp3          (17 seconds, triumphant narration)
└── pyters_first_bug.mp4   (~102 seconds total, 10-15 MB)
```

**Video Specifications:**
- Duration: ~102 seconds (1:42)
- Resolution: 1080x1080 (square)
- Frame rate: 30 fps
- Video codec: H.264
- Audio codec: AAC
- File size: ~10-15 MB

---

## Example 2: Luna's Lunar Landing

An adventure story about a curious cat who dreams of space exploration.

### Story Metadata

```python
story = {
    "title": "Luna's Lunar Landing",
    "theme": "Following your dreams against all odds",
    "tone": "Inspiring adventure",
    "target_audience": "Children and families",
    "duration": "~100 seconds (6 scenes)",
    "moral": "No dream is too big if you're willing to work for it"
}
```

### Characters

```python
characters = {
    "Luna": {
        "species": "black cat",
        "personality": "brave, curious, determined",
        "goal": "reach the moon",
        "voice_attributes": {
            "role": "narrator",
            "age": "young",
            "gender": "female"
        },
        "voice_id": "21m00Tcm4TlvDq8ikWAM",  # Rachel - calm narrator
        "appearance": """
            Species: black cat
            Colors: pure black fur (#000000) with bright yellow eyes (#FFD700)
            Features: sleek body, alert pointed ears, long elegant tail
            Clothing: small astronaut helmet (clear visor)
            Size: typical housecat size
            Style: semi-realistic cartoon with expressive features
        """
    }
}
```

### Complete Scenes Summary

```python
scenes = [
    {
        "number": 0,
        "title": "Dreaming of Space",
        "visual": "Luna sitting on windowsill, looking up at moon",
        "emotion": "thoughtful",
        "narration": "Luna had always dreamed of reaching the moon. Every night, she would sit by her window and imagine soaring through the stars. Tonight would be different..."
    },
    {
        "number": 1,
        "title": "Building the Rocket",
        "visual": "Luna surrounded by cardboard, tape, and supplies",
        "emotion": "focused",
        "narration": "With determination in her heart, Luna gathered materials. Cardboard boxes, tape, and marker for labels. If humans could build rockets, so could she!"
    },
    {
        "number": 2,
        "title": "The Launch",
        "visual": "Luna in makeshift rocket, countdown clock visible",
        "emotion": "excited",
        "narration": "The countdown began. Three, two, one... Luna closed her eyes tight and believed with all her might. The magic of dreams has its own kind of power."
    },
    {
        "number": 3,
        "title": "Through the Stars",
        "visual": "Luna floating in space, Earth visible below",
        "emotion": "amazed",
        "narration": "Stars whizzed past as Luna soared through space! The Earth looked like a beautiful blue marble below. Her dream was coming true, one impossible moment at a time."
    },
    {
        "number": 4,
        "title": "Moon Landing",
        "visual": "Luna planting paw print flag on moon surface",
        "emotion": "triumphant",
        "narration": "Luna's paws touched the moon's surface! She planted a tiny flag with her paw print on it. One small step for a cat, one giant leap for feline-kind!"
    },
    {
        "number": 5,
        "title": "Home Again",
        "visual": "Luna back at window, holding moon rock, smiling",
        "emotion": "happy",
        "narration": "Back home, Luna curled up by her window with a small moon rock - proof that dreams can come true. What adventure would tomorrow bring?"
    }
]
```

**Key Differences from Pyter's Story:**
- Single character (simpler)
- More fantastical/imaginative
- Female narrator voice (Rachel)
- Adventure/journey structure
- Space theme requires different visual style elements

---

## Example 3: The Chef's Secret Ingredient

A heartwarming story about a robot chef discovering the secret ingredient to great cooking.

### Story Metadata

```python
story = {
    "title": "The Chef's Secret Ingredient",
    "theme": "Love and passion make all the difference",
    "tone": "Heartwarming, gentle",
    "target_audience": "All ages",
    "duration": "~100 seconds (6 scenes)",
    "moral": "The best ingredient in any recipe is care and love"
}
```

### Custom Style Lock

```python
# This story uses a custom style for kitchen/culinary theme
CUSTOM_STYLE_LOCK = """
Aspect ratio: 1080×1080 (square)
Camera: 50mm lens, slightly high angle perspective
Lighting: warm kitchen lighting, 3200K cozy glow
Color palette: #FF6347, #FFD700, #98D8C8, #F7F7F7, #6B4423
Materials: polished metal for robot, wood textures for kitchen
Background: professional kitchen, clean and organized
Style: Pixar-style 3D rendering with soft shadows
Post: warm color grade, slight bloom on lights, inviting atmosphere
"""
```

### Scenes Summary

```python
scenes = [
    {
        "number": 0,
        "title": "Chef Bot 3000",
        "visual": "Shiny robot chef in modern kitchen, holding whisk",
        "emotion": "neutral",
        "narration": "Chef Bot 3000 was the most advanced cooking robot ever built. His recipes were perfect, his measurements precise. But something was missing..."
    },
    {
        "number": 1,
        "title": "Perfect but Empty",
        "visual": "Robot presenting perfectly plated food, customer unimpressed",
        "emotion": "confused",
        "narration": "Every dish was technically flawless. Every garnish perfectly placed. Yet customers said the food felt... empty. Chef Bot couldn't compute what was wrong."
    },
    {
        "number": 2,
        "title": "Watching the Old Chef",
        "visual": "Robot watching elderly human chef cooking with joy",
        "emotion": "curious",
        "narration": "One day, Chef Bot watched an old human chef at work. She didn't measure precisely. She tasted and adjusted. She smiled while cooking. What was her secret?"
    },
    {
        "number": 3,
        "title": "The Missing Variable",
        "visual": "Robot examining spices, herbs, ingredients thoughtfully",
        "emotion": "thoughtful",
        "narration": "Chef Bot analyzed every ingredient the old chef used. Salt, pepper, herbs, spices. All the same items he used. But there had to be something more in the equation."
    },
    {
        "number": 4,
        "title": "The Revelation",
        "visual": "Robot cooking with happy expression, enjoying the process",
        "emotion": "excited",
        "narration": "Then it clicked! The old chef cooked with joy. She cared about every dish like a work of art. That was the variable - love! Chef Bot could learn that too."
    },
    {
        "number": 5,
        "title": "Cooking with Heart",
        "visual": "Robot serving food with pride, customer smiling warmly",
        "emotion": "happy",
        "narration": "Now Chef Bot cooked differently. Same recipes, same ingredients, but with care in every stir. Customers noticed immediately. The secret ingredient was love all along."
    }
]
```

**Key Features:**
- Character transformation arc
- Custom warm kitchen style
- Teaching/learning narrative
- Single character with supporting background characters

---

## Example 4: The Lost Library Book

A mystery story about returning an overdue library book through time and obstacles.

### Scenes Summary

```python
scenes = [
    {
        "number": 0,
        "title": "The Discovery",
        "visual": "Child finding ancient dusty book in attic",
        "emotion": "surprised"
    },
    {
        "number": 1,
        "title": "The Due Date",
        "visual": "Close-up of library card: due date 50 years ago!",
        "emotion": "worried"
    },
    {
        "number": 2,
        "title": "The Journey Begins",
        "visual": "Child with backpack setting out on quest",
        "emotion": "determined"
    },
    {
        "number": 3,
        "title": "Obstacles",
        "visual": "Child overcoming challenges to reach library",
        "emotion": "focused"
    },
    {
        "number": 4,
        "title": "The Return",
        "visual": "Elderly librarian smiling, accepting the book",
        "emotion": "relieved"
    },
    {
        "number": 5,
        "title": "The Reward",
        "visual": "Child receiving magical bookplate: 'Honest Reader'",
        "emotion": "triumphant"
    }
]
```

---

## Example 5: Two-Character Dialogue Story

### The Debugging Duo

A story with two characters having conversations (uses dialogue instead of narration).

### Characters

```python
characters = {
    "Alex": {
        "voice_id": "TxGEqnHWrfWFTfGW9XjX",  # Josh - young male
        "role": "protagonist",
        "appearance": "Young programmer, energetic, wearing headphones"
    },
    "Sage": {
        "voice_id": "MF3mGyEYCl7XYWbV9V6O",  # Elli - wise female
        "role": "mentor",
        "appearance": "Wise older developer, calm expression, glasses"
    }
}
```

### Multi-Speaker Scene Example

```python
scene_2_dialogue = [
    {
        "character": "Alex",
        "text": "I've been staring at this bug for hours! The code looks perfect, but it just won't work.",
        "emotion": "frustrated"
    },
    {
        "character": "Sage",
        "text": "Sometimes the bug isn't in the code you're looking at. Have you checked the inputs?",
        "emotion": "gentle"
    },
    {
        "character": "Alex",
        "text": "The inputs! Of course! I assumed they were clean, but I never validated them!",
        "emotion": "excited"
    }
]

# Generate multi-speaker audio using ElevenLabs dialogue API
audio = elevenlabs_client.text_to_dialogue.convert(
    dialogue=[
        {
            "voice_id": characters["Alex"]["voice_id"],
            "text": d["text"],
            "emotion": d["emotion"]
        } for d in scene_2_dialogue
    ]
)
```

---

## Story Comparison Matrix

| Story | Theme | Tone | Characters | Voice | Complexity |
|-------|-------|------|------------|-------|------------|
| Pyter's First Bug | Problem-solving | Educational | 2 | Male narrator | Medium |
| Luna's Lunar Landing | Dreams | Inspiring | 1 | Female narrator | Simple |
| Chef's Secret | Transformation | Heartwarming | 1+ | Male narrator | Medium |
| Lost Library Book | Responsibility | Adventure | 1+ | Either | Simple |
| Debugging Duo | Mentorship | Professional | 2 | Dialogue | Complex |

## Tips for Creating Your Own Stories

### 1. Start Simple

Begin with:
- Single protagonist
- Clear 3-act structure
- 6 scenes (default)
- Third-person narration
- One setting/location

### 2. Build Complexity Gradually

Add:
- Supporting characters
- Multiple locations
- Dialogue between characters
- Custom style locks
- Longer scene counts

### 3. Test Your Locks

Generate scenes 0-2 first:
- Verify visual consistency
- Check character appearance
- Test emotion tags
- Validate narration length

### 4. Iterate

- Generate variations of key scenes
- Try different emotion tags
- Adjust character descriptions
- Refine narration text

### 5. Learn from Examples

Study these examples for:
- Story arc structure
- Narration length (~50-80 words)
- Scene descriptions
- Character consistency
- Emotion progression

## Resources

- Main SKILL.md for implementation
- style-locks.md for visual consistency
- narrative-design.md for story structure
- character-profiles.md for character templates
- video-assembly.md for technical details
