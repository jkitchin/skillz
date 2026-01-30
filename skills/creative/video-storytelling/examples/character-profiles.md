# Character Profile Templates

## Overview

Character profiles define the visual appearance, personality, and voice characteristics for story characters. Consistent character definitions are essential for maintaining visual coherence across multi-scene video stories.

## Character Profile Template

```python
character = {
    # Basic Information
    "name": "CharacterName",
    "species": "What type of character (human, animal, robot, etc.)",
    "role": "protagonist | supporting | antagonist | narrator",

    # Personality
    "personality": "Key traits (curious, brave, thoughtful, etc.)",
    "goal": "What the character wants",
    "arc": "How the character changes (optional)",

    # Voice Attributes (for voice assignment)
    "voice_attributes": {
        "role": "narrator | character",
        "age": "young | middle | old",
        "gender": "male | female | neutral"
    },

    # ElevenLabs Voice
    "voice_id": "ElevenLabs voice ID",
    "voice_name": "Voice name for reference",

    # Visual Appearance (CHARACTER_LOCK)
    "appearance": """
        Species: [detailed species/type description]
        Colors: [primary colors with hex codes]
        Features: [distinctive visual characteristics]
        Clothing: [what they wear, if anything]
        Size: [relative size and scale]
        Age: [visual age appearance]
        Expression: [default expression or mood]
        Style: [art style, should match global STYLE_LOCK]
    """
}
```

## Example Character Profiles

### Example 1: Protagonist - Adventurous Child

```python
character_profiles = {
    "Mira": {
        "name": "Mira",
        "species": "human child",
        "role": "protagonist",

        "personality": "curious, brave, imaginative, kind-hearted",
        "goal": "explore the magical forest and find the lost star",
        "arc": "learns that courage means being afraid but doing it anyway",

        "voice_attributes": {
            "role": "character",
            "age": "young",
            "gender": "female"
        },

        "voice_id": "EXAVITQu4vr4xnSDxMaL",  # Bella
        "voice_name": "Bella (young expressive female)",

        "appearance": """
            Species: human child, approximately 8 years old
            Colors: brown skin tone, dark curly hair (#2C1810),
                    bright red jacket (#E74C3C), blue jeans (#3498DB)
            Features: large expressive brown eyes, friendly smile,
                      curly hair in two puffs, small backpack
            Clothing: red adventure jacket with pockets, comfortable jeans,
                      hiking boots, small yellow backpack
            Size: child-sized, about 4 feet tall
            Age: 8 years old in appearance
            Expression: curious and determined, slight smile
            Style: semi-realistic cartoon with gentle shading and clear features
        """
    }
}
```

### Example 2: Supporting Character - Wise Mentor

```python
character_profiles["Professor Oak"] = {
    "name": "Professor Oak",
    "species": "elderly human scientist",
    "role": "supporting",

    "personality": "wise, patient, encouraging, knowledgeable",
    "goal": "help young scientists discover their potential",

    "voice_attributes": {
        "role": "character",
        "age": "old",
        "gender": "male"
    },

    "voice_id": "ErXwobaYiN019PkySvjV",  # Antoni
    "voice_name": "Antoni (older wise male)",

    "appearance": """
        Species: human male, elderly scientist
        Colors: white hair and beard (#F0F0F0), tan lab coat (#D4A574),
                brown vest (#8B4513), grey pants (#708090)
        Features: kind wrinkled face, round spectacles,
                  warm smile, slightly hunched posture
        Clothing: traditional tan lab coat, brown vest underneath,
                  bow tie, pocket protector with pens
        Size: average height adult male, slightly shorter with age
        Age: elderly, approximately 70 years old
        Expression: warm grandfatherly smile, twinkling eyes
        Style: semi-realistic cartoon with expressive wrinkles and soft features
    """
}
```

### Example 3: Animal Character - Loyal Companion

```python
character_profiles["Rufus"] = {
    "name": "Rufus",
    "species": "golden retriever dog",
    "role": "supporting",

    "personality": "loyal, playful, protective, enthusiastic",
    "goal": "stay by his owner's side through all adventures",

    "voice_attributes": {
        "role": "narrator",  # Dogs don't talk, but narration can express thoughts
        "age": "young",
        "gender": "male"
    },

    "voice_id": "TxGEqnHWrfWFTfGW9XjX",  # Josh
    "voice_name": "Josh (energetic young male for thoughts)",

    "appearance": """
        Species: golden retriever dog
        Colors: golden tan fur (#D4A017), white chest patch (#FFFFFF),
                brown nose (#8B4513), dark brown eyes (#3E2723)
        Features: floppy ears, wagging tail, friendly panting expression,
                  shiny coat, alert posture
        Clothing: red collar with silver tag
        Size: large dog, comes up to adult's waist
        Age: young adult dog, energetic and fit
        Expression: happy panting smile, bright attentive eyes
        Style: semi-realistic cartoon with soft fur texture
    """
}
```

### Example 4: Anthropomorphic Character - Robot

```python
character_profiles["Byte"] = {
    "name": "Byte",
    "species": "small helpful robot",
    "role": "supporting",

    "personality": "helpful, logical, learning emotions, cheerful",
    "goal": "assist humans and understand feelings",

    "voice_attributes": {
        "role": "character",
        "age": "young",
        "gender": "neutral"
    },

    "voice_id": "pNInz6obpgDQGcFmaJgB",  # Adam
    "voice_name": "Adam (neutral friendly voice)",

    "appearance": """
        Species: small humanoid robot
        Colors: white chassis (#FFFFFF), blue accent panels (#3498DB),
                glowing cyan eyes (#00FFFF), silver joints (#C0C0C0)
        Features: rounded head with digital display face,
                  articulated arms and legs, antenna on head,
                  friendly glowing eyes that show emotions
        Clothing: built-in appearance, no separate clothing
        Size: small, about 2 feet tall (child-sized)
        Age: appears new/modern in design
        Expression: digital smile on display face, bright friendly eyes
        Style: clean 3D render style with slight cartoon proportions
    """
}
```

### Example 5: Fantasy Character - Magical Creature

```python
character_profiles["Sparkle"] = {
    "name": "Sparkle",
    "species": "fairy dragon",
    "role": "supporting",

    "personality": "mischievous, magical, protective, whimsical",
    "goal": "protect the enchanted forest",

    "voice_attributes": {
        "role": "character",
        "age": "young",
        "gender": "female"
    },

    "voice_id": "EXAVITQu4vr4xnSDxMaL",  # Bella
    "voice_name": "Bella (bright playful female)",

    "appearance": """
        Species: small fairy dragon (magical creature)
        Colors: iridescent purple scales (#9B59B6),
                pink wing membranes (#FF69B4),
                golden horns and claws (#FFD700),
                sparkle effects (#FFFFFF with glow)
        Features: large expressive eyes, delicate translucent wings,
                  small horns, long graceful tail with tuft,
                  magical sparkles around body
        Clothing: natural appearance, no clothing
        Size: very small, fits in human hand (6 inches long)
        Age: young but ancient species
        Expression: playful smirk, twinkling eyes, wings slightly spread
        Style: fantasy illustration with magical glow effects
    """
}
```

## Voice Assignment Reference

### ElevenLabs Voice Catalog

**Narrative Voices (Third-Person):**

```python
NARRATOR_VOICES = {
    "George": {
        "voice_id": "JBFqnCBsd6RMkjVDRZzb",
        "gender": "male",
        "age": "adult",
        "description": "Warm professional narrator, perfect for audiobooks"
    },
    "Rachel": {
        "voice_id": "21m00Tcm4TlvDq8ikWAM",
        "gender": "female",
        "age": "young adult",
        "description": "Calm soothing narrator, great for gentle stories"
    }
}
```

**Character Voices (Dialogue):**

```python
CHARACTER_VOICES = {
    # Young Characters
    "Josh": {
        "voice_id": "TxGEqnHWrfWFTfGW9XjX",
        "gender": "male",
        "age": "young",
        "description": "Energetic upbeat young male"
    },
    "Bella": {
        "voice_id": "EXAVITQu4vr4xnSDxMaL",
        "gender": "female",
        "age": "young",
        "description": "Bright expressive young female"
    },

    # Adult Characters
    "Arnold": {
        "voice_id": "VR6AewLTigWG4xSOukaG",
        "gender": "male",
        "age": "middle",
        "description": "Confident middle-aged male"
    },
    "Elli": {
        "voice_id": "MF3mGyEYCl7XYWbV9V6O",
        "gender": "female",
        "age": "middle",
        "description": "Warm middle-aged female"
    },

    # Elderly Characters
    "Antoni": {
        "voice_id": "ErXwobaYiN019PkySvjV",
        "gender": "male",
        "age": "old",
        "description": "Wise older male, great for mentors"
    },
    "Dorothy": {
        "voice_id": "ThT5KcBeYPX3keUQqHPh",
        "gender": "female",
        "age": "old",
        "description": "Kind older female, grandmother voice"
    },

    # Neutral/Flexible
    "Adam": {
        "voice_id": "pNInz6obpgDQGcFmaJgB",
        "gender": "neutral",
        "age": "young-adult",
        "description": "Friendly neutral voice, good for robots/AI"
    }
}
```

## Character Consistency Checklist

When creating character profiles, ensure:

- [ ] **Clear species definition** - Specific type, not generic
- [ ] **Exact colors with hex codes** - Prevents color drift
- [ ] **Distinctive features listed** - Makes character recognizable
- [ ] **Consistent size/scale** - Maintains proportions across scenes
- [ ] **Appropriate voice match** - Voice fits character's age/personality
- [ ] **Style matches global STYLE_LOCK** - Visual coherence
- [ ] **Detailed enough for regeneration** - Someone else could recreate
- [ ] **Default expression defined** - Starting emotional state

## Multi-Character Story Example

### The Magical Garden Adventure

```python
story_characters = {
    "Emma": {
        "name": "Emma",
        "species": "human child",
        "role": "protagonist",
        "personality": "curious, kind, gentle with plants",
        "voice_id": "EXAVITQu4vr4xnSDxMaL",  # Bella
        "appearance": """
            Species: human child, 7 years old
            Colors: light skin, blonde hair (#F4D03F), green dress (#27AE60)
            Features: braided hair with flower clips, freckles, bright smile
            Clothing: green sundress, white sandals, sun hat
            Size: small child, approximately 3.5 feet tall
            Expression: gentle curious smile
            Style: soft cartoon with warm lighting
        """
    },

    "Sprout": {
        "name": "Sprout",
        "species": "sentient plant creature",
        "role": "supporting",
        "personality": "shy, wise, speaks in plant metaphors",
        "voice_id": "ThT5KcBeYPX3keUQqHPh",  # Dorothy
        "appearance": """
            Species: magical plant creature, living flower
            Colors: green stem body (#2ECC71), pink flower petals (#FF69B4),
                    yellow center (#F1C40F)
            Features: flower head as face, vine arms, root legs,
                      expressive petal arrangement
            Clothing: natural plant form
            Size: small, Emma's knee height
            Expression: shy gentle smile
            Style: whimsical cartoon with soft glow
        """
    },

    "Digger": {
        "name": "Digger",
        "species": "friendly garden mole",
        "role": "supporting",
        "personality": "hardworking, practical, underground expert",
        "voice_id": "TxGEqnHWrfWFTfGW9XjX",  # Josh
        "appearance": """
            Species: cartoon mole
            Colors: dark brown fur (#5D4037), pink nose (#FFB6C1),
                    small black eyes (#000000)
            Features: large digging claws, soft fuzzy fur,
                      tiny glasses, dirt smudges
            Clothing: small yellow hard hat, tool belt
            Size: very small, fits in Emma's hands
            Expression: cheerful determined look
            Style: cute cartoon with texture details
        """
    }
}
```

## Character Interaction Guidelines

### Visual Positioning

When multiple characters appear in a scene:

```python
scene_visual = """
Emma kneels in the garden (left side of frame),
Sprout stands beside her at ground level (center),
Digger peeks out from a small hole (right foreground).
All three are facing slightly toward each other,
creating a triangular composition. Maintain all
character appearances as defined in CHARACTER_LOCK.
"""
```

### Voice Dialogue

```python
scene_dialogue = [
    {
        "character": "Emma",
        "voice_id": story_characters["Emma"]["voice_id"],
        "text": "Look! The flowers are wilting. Can you help, Sprout?",
        "emotion": "concerned"
    },
    {
        "character": "Sprout",
        "voice_id": story_characters["Sprout"]["voice_id"],
        "text": "They need water, dear child. Their roots thirst for rain.",
        "emotion": "gentle"
    },
    {
        "character": "Digger",
        "voice_id": story_characters["Digger"]["voice_id"],
        "text": "I can dig channels to bring water! Let's work together!",
        "emotion": "excited"
    }
]
```

## Character Evolution (Optional)

Characters can change appearance over the story (with intention):

```python
# Initial appearance
character_start = """
    Clothing: clean white shirt, neat hair
    Expression: uncertain, nervous
"""

# After transformation scene
character_end = """
    Clothing: same white shirt but with dirt and grass stains
              (shows adventure happened)
    Expression: confident smile, proud posture
"""
```

**Important:** Only change appearance intentionally as part of story progression.

## Tips for Great Characters

### 1. Make Them Distinctive

**Bad:**
```
A cat with fur and eyes
```

**Good:**
```
A sleek black cat with bright yellow eyes (#FFD700),
white chest patch, pink nose, always alert ears
```

### 2. Match Voice to Character

- Young energetic character → Josh or Bella
- Wise mentor → Antoni or Dorothy
- Neutral narrator → George or Rachel
- Cute creature → Bella with higher pitch

### 3. Define Personality Through Appearance

- **Curious character:** Wide eyes, forward lean, open expression
- **Wise character:** Gentle smile, calm posture, thoughtful gaze
- **Brave character:** Strong stance, determined expression, ready posture
- **Shy character:** Smaller posture, softer features, hesitant look

### 4. Consistency is Key

Use exact same description in every scene's CHARACTER_LOCK section.

### 5. Test Early

Generate scenes 0 and 1 to verify:
- Character looks consistent
- Colors are correct
- Features are recognizable
- Voice sounds appropriate

## Resources

- Main SKILL.md for implementation
- style-locks.md for visual consistency
- narrative-design.md for character arcs
- example-stories.md for complete character examples
- ElevenLabs voice library: https://elevenlabs.io/voice-library
